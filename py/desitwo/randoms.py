"""
desitwo.randoms
==================

Monte Carlo imaging at the pixel level. Ported from desitarget.randoms.

.. _`Legacy Surveys bitmasks`: https://www.legacysurvey.org/dr10/bitmasks/

"""
import os
import photutils
import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from time import time

from desitarget.randoms import randoms_in_a_brick_from_edges, bricklookup, get_dust
from desitarget.skyfibers import get_brick_info

# ADM the parallelization script.
from desitarget.internal import sharedmem

# ADM set up the DESI default logger.
from desiutil.log import get_logger
log = get_logger()

# ADM start the clock.
start = time()


def is_in_gaia_mask_arjun(ras, decs, clip=True,
                          fn='/global/homes/a/arjundey/ODIN/bright_stars_for_mask.fits'):
    """Mask bright stars using Arjun's Gaia mask. Original version.

    Parameters
    ----------
    ras : :class:`~numpy.array`
        Right Ascensions of interest (degrees).
    decs : :class:`~numpy.array`
        Declinations of interest (degrees).
    clip : :class:`bool`
        If ``True`` then clip the centers of the bright star masks to be
        limited to the boundary box formed by the input `ras` and `decs`.
        Arjun originally clipped (perhaps as a speed-up).
    fn : :class:`str`
        Filename of bright star mask.

    Returns
    -------
    :class:`~numpy.array`
        An array of Booleans that is ``True`` for locations in the mask
        and ``False`` for locations outside the mask.
    :class:`str`
        Name of the filename used/passed to read the mask.

    Notes
    -----
    - This is a close copy of Arjun Dey's code from a notebook at
      /global/homes/a/arjundey/ODIN/ODIN_COSMOS_N419.ipynb
    - :func:`is_in_gaia_mask()` should produce identical results, and
      should be faster.
    """
    start = time()

    # Read in the Bright star mask
    # Bright stars with G < 9 mag from Gaia DR2.
    br = Table.read(fn, hdu=1)

    # ADM populate needed info in one place.
    dt = [('ra', '>f8'), ('dec', 'f8'), ('brightstarflag', '<i4')]
    alldatals = np.zeros(len(ras), dtype=dt)
    alldatals['ra'] = ras
    alldatals['dec'] = decs

    # Restrict bright star mask to size of the ODIN catalog and to.
    if clip:
        i = ((br['RA'] >= np.min(alldatals['ra']))
             & (br['RA'] <= np.max(alldatals['ra']))
             & (br['DEC'] >= np.min(alldatals['dec']))
             & (br['DEC'] <= np.max(alldatals['dec'])))
        br = br[i]

    # ADM also clip to appropriate G-magnitude.
    ii = br['GAIA_PHOT_G_MEAN_MAG'] <= 9.0
    br = br[ii]

    # Match the ODIN LAE candidates with the bright star catalog.
    c1 = SkyCoord(ra=alldatals['ra']*u.degree, dec=alldatals['dec']*u.degree)
    c2 = SkyCoord(ra=br['RA']*u.degree, dec=br['DEC']*u.degree)
    idx, d2d, d3d = c1.match_to_catalog_sky(c2)
    max_sep = 1.0*u.arcsec
    max_sep = 0.07*(6.3/br['GAIA_PHOT_G_MEAN_MAG'])**2

    # ADM I changed this print statement to allow logging.
    log.info(f'Number of bright stars = {len(br)}...t={time()-start:.1f}s')

    # Flag the sources that fall near bright stars.
    lenbr = len(br)
    badtot = 0
    # For each star in the bright star list flag the catalog sources
    # which lie within the appropriate radius.
    for i in range(lenbr):
        c2 = SkyCoord(ra=br['RA'][i]*u.degree, dec=br['DEC'][i]*u.degree)
        # matches to the catalog.
        # Old formula: msep_i = 0.07*(6.3/br['GAIA_PHOT_G_MEAN_MAG'][i])**2
        # Variable matching criterion based on the number
        msep_i = 0.07*(6.3/br['GAIA_PHOT_G_MEAN_MAG'][i])**2
        d2d = c1.separation(c2)
        ibad_i = d2d < msep_i*u.deg
        log.info('Mag= %f ; Separation= %f; Number of flagged sources = %i'
                 % (br['GAIA_PHOT_G_MEAN_MAG'][i], msep_i, len(d2d[ibad_i])))
        badtot = badtot + len(d2d[ibad_i])
        if len(d2d[ibad_i]) > 0:
            alldatals['brightstarflag'][ibad_i] = 1

    ibad = (alldatals['brightstarflag'] > 0)
    log.info(f"Total number of catalog sources flagged due to proximity "
             f"to bright stars = {np.sum(ibad)}...t={time()-start:.1f}s")

    return ibad, fn


def is_in_gaia_mask(ras, decs, clip=True,
                    fn='/global/homes/a/arjundey/ODIN/bright_stars_for_mask.fits'):
    """Mask bright stars using Arjun's Gaia mask. Quick/compact version.

    Parameters
    ----------
    ras : :class:`~numpy.array`
        Right Ascensions of interest (degrees).
    decs : :class:`~numpy.array`
        Declinations of interest (degrees).
    clip : :class:`bool`
        If ``True`` then clip the centers of the bright star masks to be
        limited to the boundary box formed by the input `ras` and `decs`.
        Potentially useful as a speed-up.
    fn : :class:`str`
        Filename of bright star mask.

    Returns
    -------
    :class:`~numpy.array`
        An array of Booleans that is ``True`` for locations in the mask
        and ``False`` for locations outside the mask.
    :class:`str`
        Name of the filename used/passed to read the mask.

    Notes
    -----
     - :func:`is_in_gaia_mask_arjun()` should produce identical results,
       but this code should be faster.
    """
    start = time()

    # ADM Read in the bright star mask from Arjun's home directory.
    Mx = Table.read(fn)

    # ADM Limit bright star mask to size of ODIN catalog, if requested.
    if clip:
        ii = ((Mx['RA'] >= np.min(ras)) & (Mx['RA'] <= np.max(ras))
              & (Mx['DEC'] >= np.min(decs)) & (Mx['DEC'] <= np.max(decs)))
        Mx = Mx[ii]
    # ADM Also restrict bright star mask to appropriate magnitude limit.
    ii = Mx['GAIA_PHOT_G_MEAN_MAG'] <= 9.0
    Mx = Mx[ii]

    log.info(f"Working with {len(Mx)} bright star masks t={time()-start:.1f}s")

    # ADM if there are no masks in the area, exit with zero masking.
    if len(Mx) == 0:
        log.info("So no masking is required")
        inmask = np.zeros(len(ras), dtype="?")
        return inmask, fn

    # ADM Match the passed locations with the bright star mask.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
    cMx = SkyCoord(ra=Mx['RA']*u.degree, dec=Mx['DEC']*u.degree)
    # ADM Organize the sense of the match so each passed location has an
    # ADM associated closest bright star mask (identified by iiMx)..
    iiMx, sep, _ = c.match_to_catalog_sky(cMx)

    # ADM The maximum separation for each mask.
    maxsep = 0.07*(6.3/Mx['GAIA_PHOT_G_MEAN_MAG'])**2

    # ADM A coordinate is in a mask if the separation from the mask is
    # ADM less than the maximum separation associated with the mask.
    inmask = sep < maxsep[iiMx]*u.deg

    log.info(f'{np.sum(inmask)} locations masked by stars...t={time()-start:.1f}s')

    return inmask, fn


def quantities_at_positions_in_a_brick(ras, decs, brickname, direc,
                                       aprad=0.75):
    """Observational quantities (per-band) at positions in a Legacy Surveys brick.

    Parameters
    ----------
    ras : :class:`~numpy.array`
        Right Ascensions of interest (degrees).
    decs : :class:`~numpy.array`
        Declinations of interest (degrees).
    brickname : :class:`str`
        Name of brick which contains RA/Dec positions, e.g., '1351p320'.
    direc : :class:`str`
       The root directory pointing to a set of Legacy-Surveys-like imaging
       files, e.g., /global/cfs/cdirs/cosmo/work/users/dstn/ODIN/xmm-N419/coadd/
    aprad : :class:`float`, optional, defaults to 0.75
        Radii in arcsec of aperture for which to derive sky/fiber fluxes.
        Defaults to the DESI fiber radius. If aprad < 1e-8 is passed,
        the code to produce these values is skipped, as a speed-up, and
        `apflux_` output values are set to zero.

    Returns
    -------
    :class:`dictionary`
       The number of observations (`nobsx`), PSF depth (`psfdepth_x`)
       galaxy depth (`galdepth_x`), PSF size (`psfsize_x`), sky
       background (`apflux_x`) and inverse variance (`apflux_ivar_x`)
       at each passed position in each band x. Plus, the `maskbits`
       information at each passed position for the brick.
       Also adds a unique `objid` for each random.
    """
    # ADM we're always reading files relative to extension number 1.
    extn_nb = 1

    # ADM guard against too low a density of random locations.
    npts = len(ras)
    if npts == 0:
        msg = 'brick {} is empty. Increase the density of random points!'.format(brickname)
        log.critical(msg)
        raise ValueError(msg)

    # ADM a list to populate with the files required to run the code.
    fnlist = []

    # ADM the output dictionary.
    qdict = {}

    # as a speed up, we assume all images in different filters for the brick have the same WCS
    # -> if we have read it once (iswcs=True), we use this info
    iswcs = False
    # ADM this will store the instrument name the first time we touch the wcs
    instrum = None

    # ADM Some choices for filter names. Default to ODIN filters...
    filters = ['N419', 'N501', 'N673']
    # ADM ...different filters if we're using Suprime-Cam.
    if "suprime" in direc:
        filters = ['I-A-L427', 'I-A-L464', 'I-A-L484', 'I-A-L505', 'I-A-L527']

    rootdir = os.path.join(direc, 'coadd', brickname[:3], brickname)
    fileform = os.path.join(rootdir, 'legacysurvey-{}-{}-{}.fits.fz')
    # ADM loop through the filters and store the number of observations
    # ADM etc. at the RA and Dec positions of the passed points.
    for filt in filters:
        # ADM the input file labels, and output column names and output
        # ADM formats for each of the quantities of interest.
        qnames = zip(['nexp', 'depth', 'galdepth', 'psfsize', 'image'],
                     ['nobs', 'psfdepth', 'galdepth', 'psfsize', 'apflux'],
                     ['i2', 'f4', 'f4', 'f4', 'f4'])
        for qin, qout, qform in qnames:
            fn = fileform.format(brickname, qin, filt)
            # ADM only process the WCS if there's a file for this filter.
            # ADM also skip calculating aperture fluxes if aprad ~ 0.
            if os.path.exists(fn) and not (qout == 'apflux' and aprad < 1e-8):
                img = fits.open(fn)[extn_nb]
                if not iswcs:
                    w = WCS(img.header)
                    x, y = w.all_world2pix(ras, decs, 0)
                    iswcs = True
                # ADM get the quantity of interest at each location and
                # ADM store in a dictionary with the filter and quantity.
                if qout == 'apflux':
                    # ADM special treatment to photometer sky.
                    # ADM Read in the ivar image.
                    fnivar = fileform.format(brickname, 'invvar', filt)
                    ivar = fits.open(fnivar)[extn_nb].data
                    with np.errstate(divide='ignore', invalid='ignore'):
                        # ADM ivars->errors, guard against 1/0.
                        imsigma = 1./np.sqrt(ivar)
                        imsigma[ivar == 0] = 0
                    # ADM aperture photometry at requested radius (aprad).
                    apxy = np.vstack((x, y)).T
                    aper = photutils.CircularAperture(apxy, aprad)
                    p = photutils.aperture_photometry(img.data, aper, error=imsigma)
                    # ADM store the results.
                    qdict[qout+'_'+filt] = np.array(p.field('aperture_sum'))
                    err = p.field('aperture_sum_err')
                    with np.errstate(divide='ignore', invalid='ignore'):
                        # ADM errors->ivars, guard against 1/0.
                        ivar = 1./err**2.
                        ivar[err == 0] = 0.
                    qdict[qout+'_ivar_'+filt] = np.array(ivar)
                else:
                    qdict[qout+'_'+filt] = img.data[y.round().astype("int"), x.round().astype("int")]
            # ADM if the file doesn't exist, set quantities to zero.
            else:
                if qout == 'apflux':
                    qdict['apflux_ivar_'+filt] = np.zeros(npts, dtype=qform)
                qdict[qout+'_'+filt] = np.zeros(npts, dtype=qform)

    # ADM add the MASKBITS information.
    fn = os.path.join(rootdir,
                      'legacysurvey-{}-maskbits.fits.fz'.format(brickname))
    # ADM only process the WCS if there's a file for this filter.
    mnames = zip([extn_nb],
                 ['maskbits'],
                 ['>i2'])
    for mextn, mout, mform in mnames:
        if os.path.exists(fn):
            img = fits.open(fn)[mextn]
            # ADM use the WCS for the per-filter quantities if it exists.
            if not iswcs:
                w = WCS(img.header)
                x, y = w.all_world2pix(ras, decs, 0)
                iswcs = True
            # ADM add the maskbits to the dictionary.
            qdict[mout] = img.data[y.round().astype("int"), x.round().astype("int")]
        else:
            # ADM if no files are found, populate with zeros.
            qdict[mout] = np.zeros(npts, dtype=mform)
            # ADM if there was no maskbits file, populate with BAILOUT.
            if mout == 'maskbits':
                qdict[mout] |= 2**10

    # ADM assign a reasonable OBJID based on order in RA.
    qdict["objid"] = np.argsort(ras)

    return qdict


def get_quantities_in_a_brick(ramin, ramax, decmin, decmax, brickname, direc,
                              density=100000, dustdir=None, aprad=0.75, seed=1):
    """NOBS, DEPTHS etc. (per-band) for random points in a brick.

    Parameters
    ----------
    ramin, ramax, decmin, decmax : :class:`float`
        Minimum and maximum brick "edges" in RA and Dec.
    brickname : :class:`~numpy.array`
        Brick name that corresponds to the brick edges, e.g., '1351p320'.
    direc : :class:`str`
       The root directory of a set of Legacy-Surveys-like imaging files
       e.g., /global/cfs/cdirs/cosmo/work/users/dstn/ODIN/xmm-N419/coadd/
    density : :class:`int`, optional, defaults to 100,000
        Number of random points to return per sq. deg. As a typical brick
        is ~0.25 x 0.25 sq. deg. ~0.0625*density points will be returned.
    dustdir : :class:`str`, optional, defaults to $DUST_DIR+'/maps'
        The root directory pointing to SFD dust maps. If not
        sent the code will try to use $DUST_DIR+'/maps' before failing.
    aprad : :class:`float`, optional, defaults to 0.75
        Radii in arcsec of aperture for which to derive sky/fiber fluxes.
        Defaults to the DESI fiber radius. If aprad < 1e-8 is passed,
        values are skipped as a speed-up and returned as zero.
    seed : :class:`int`, optional, defaults to 1
        Random seed used to generate RA/Dec positions in the function
        :func:`~desitarget.randoms.randoms_in_a_brick_from_edges()`.

    Returns
    -------
    :class:`~numpy.ndarray`
        OBJID:
            A unique (to each brick) source identifier.
        BRICKID, BRICKNAME:
            ID and name that correspond to the passed brick name.
        RA, DEC:
            Right Ascension, Declination of a random location.
        NOBS_X:
            Number of observations in band X.
        PSFDEPTH_X, GALDEPTH_X:
            PSF and galaxy depth at this location in band X.
        PSFSIZE_X:
            Weighted average PSF FWHM (arcsec) in band X.
        APFLUX_X:
            Sky background extracted in `aprad`. Zero for `aprad` < 1e-8.
        APFLUX_IVAR_X:
            Inverse variance of sky background. Zero for `aprad` < 1e-8.
        MASKBITS:
            Mask information. See, e.g. the `Legacy Surveys bitmasks`_.
        EBV:
            E(B-V) at this location from the SFD dust maps.
    """
    # ADM only intended to work on one brick, so die for larger arrays.
    if not isinstance(brickname, str):
        log.fatal("Only one brick can be passed at a time!")
        raise ValueError

    # ADM generate random points in the brick at the requested density.
    ras, decs = randoms_in_a_brick_from_edges(ramin, ramax, decmin, decmax,
                                              density=density, wrap=False,
                                              seed=seed)

    qdict = quantities_at_positions_in_a_brick(ras, decs, brickname,
                                               direc, aprad=aprad)

    # ADM Some choices for filter names. Default to ODIN filters...
    filters = ['N419']
    if np.any(['N501' in k or 'N673' in k for k in qdict.keys()]):
        filters += ['N501', 'N673']
    # ADM ...different filters if we're using Suprime-Cam.
    if "suprime" in direc:
        filters = ['I-A-L427', 'I-A-L464', 'I-A-L484', 'I-A-L505', 'I-A-L527']

    # ADM the dtype of the structured array to output.
    dt = [('BRICKID', '>i4'), ('BRICKNAME', 'U8'), ('OBJID', '>i4'), ('RA', '>f8'),
          ('DEC', 'f8'), ('EBV', 'f4'), ('MASKBITS', '<i4'), ('IN_ARJUN_MASK', '?')]

    for filt in filters:
        dt += [(f'NOBS_{filt}', '<i2'), (f'PSFDEPTH_{filt}', '<f4'),
               (f'GALDEPTH_{filt}', '<f4'), (f'PSFSIZE_{filt}', '<f4'),
               (f'APFLUX_{filt}', '<f8'), (f'APFLUX_IVAR_{filt}', '<f8')]

    # ADM the full structured array to output.
    qinfo = np.zeros(len(ras), dtype=dt)

    # ADM retrieve the E(B-V) values for each random point.
    ebv = get_dust(ras, decs, dustdir=dustdir)

    # ADM catch the case of a missing coadd directory.
    if len(qdict) > 0:
        # ADM store each quantity in the output structured array
        # ADM remembering the dictionary keys are lower-case.
        cols = qdict.keys()
        for col in cols:
            qinfo[col.upper()] = qdict[col]

    # ADM add the RAs/Decs, brick ID and brick name.
    brickid = bricklookup[brickname]
    qinfo["RA"], qinfo["DEC"] = ras, decs
    qinfo["BRICKNAME"], qinfo["BRICKID"] = brickname, brickid

    # ADM add the dust values.
    qinfo["EBV"] = ebv

    return qinfo


def select_randoms_bricks(brickdict, bricknames, direc, density=100000,
                          numproc=32, dustdir=None, aprad=0.75, seed=1):

    """Parallel-process a random catalog for a set of brick names.

    Parameters
    ----------
    brickdict : :class:`dict`
        Look-up dictionary for a set of bricks, as made by, e.g.
        :func:`~desitarget.skyfibers.get_brick_info()`.
    bricknames : :class:`~numpy.array`
        The names of the bricks in `brickdict` to process.
    direc : :class:`str`
       The root directory pointing to a set of Legacy-Surveys-like imaging
       files, e.g., /global/cfs/cdirs/cosmo/work/users/dstn/ODIN/xmm-N419/coadd/
    density : :class:`int`, optional, defaults to 100,000
        Number of random points to return per sq. deg. As a brick is
        ~0.25 x 0.25 sq. deg. ~0.0625*density points will be returned.
    numproc : :class:`int`, optional, defaults to 32
        The number of processes over which to parallelize.
    dustdir : :class:`str`, optional, defaults to $DUST_DIR+'maps'
        The root directory pointing to SFD dust maps. If ``None`` the
        code will try to use $DUST_DIR+'maps') before failing.
    aprad : :class:`float`, optional, defaults to 0.75
        Radii in arcsec of aperture for which to derive sky/fiber fluxes.
        Defaults to the DESI fiber radius. If aprad < 1e-8 is passed,
        values are skipped as a speed-up and returned as zero.
    seed : :class:`int`, optional, defaults to 1
        Random seed used to generate RA/Dec positions in the function
        :func:`~desitarget.randoms.randoms_in_a_brick_from_edges()`.

    Returns
    -------
    :class:`~numpy.ndarray`
        A numpy structured array with the same columns as returned by
        :func:`~desitwo.randoms.get_quantities_in_a_brick()`.
    """
    nbricks = len(bricknames)
    log.info(f"Running {nbricks} bricks from {direc} at density {density:.1e} "
             f"per sq. deg on {numproc} processors...t = {time()-start:.1f}s")

    # ADM the critical function to run on every brick.
    def _get_quantities(brickname):
        """wrapper on get_quantities_in_a_brick() given a brick name"""
        # ADM retrieve the edges for the brick that we're working on.
        bra, bdec, bramin, bramax, bdecmin, bdecmax = brickdict[brickname]

        # ADM populate the brick with random points, and retrieve the
        # ADM quantities of interest at those locations.
        randoms = get_quantities_in_a_brick(
            bramin, bramax, bdecmin, bdecmax, brickname, direc=direc,
            density=density, dustdir=dustdir, aprad=aprad, seed=seed)

        return randoms

    # ADM this is just to count bricks in _update_status.
    nbrick = np.zeros((), dtype='i8')
    t0 = time()
    # ADM write a total of 25 output messages during processing.
    interval = 1
    if nbricks >= 25:
        interval = nbricks // 25

    def _update_status(result):
        ''' wrapper function for the critical reduction operation,
            that occurs on the main parallel process '''
        if nbrick % interval == 0 and nbrick > 0:
            elapsed = time() - t0
            rate = nbrick / elapsed
            log.info('{}/{} bricks; {:.1f} bricks/sec; {:.1f} total mins elapsed'
                     .format(nbrick, nbricks, rate, elapsed/60.))
            # ADM if we're going to exceed 4 hours, warn the user.
            if nbricks/rate > 4*3600.:
                msg = 'May take > 4 hours to run. May fail on interactive nodes.'
                log.warning(msg)

        nbrick[...] += 1    # this is an in-place modification.
        return result

    # - Parallel process input files.
    if numproc > 1:
        pool = sharedmem.MapReduce(np=numproc)
        with pool:
            qinfo = pool.map(_get_quantities, bricknames, reduce=_update_status)
    else:
        qinfo = list()
        for brickname in bricknames:
            qinfo.append(_update_status(_get_quantities(brickname)))

    qinfo = np.concatenate(qinfo)

    return qinfo


def select_randoms(direc,
                   density=100000, numproc=32, dustdir=None, aprad=0.75, seed=1):
    """NOBS, DEPTHs (per-band), MASKs for random points in a Legacy Surveys DR.

    Parameters
    ----------
    direc : :class:`str`
       The root directory pointing to a set of Legacy-Surveys-like imaging
       files, e.g., /global/cfs/cdirs/cosmo/work/users/dstn/ODIN/xmm-N419/coadd/
    density : :class:`int`, optional, defaults to 100,000
        Number of random points to return per sq. deg. As a brick is
        ~0.25 x 0.25 sq. deg. ~0.0625*density points will be returned.
    numproc : :class:`int`, optional, defaults to 32
        The number of processes over which to parallelize.
    dustdir : :class:`str`, optional, defaults to $DUST_DIR+'maps'
        The root directory pointing to SFD dust maps. If None the code
        will try to use $DUST_DIR+'maps') before failing.
    aprad : :class:`float`, optional, defaults to 0.75
        Radii in arcsec of aperture for which to derive sky/fiber fluxes.
        Defaults to the DESI fiber radius. If aprad < 1e-8 is passed,
        values are skipped as a speed-up and returned as zero.
    seed : :class:`int`, optional, defaults to 1
        Random seed to use when shuffling across brick boundaries.
        The actual np.random.seed defaults to 615+`seed`. See also use
        in :func:`~desitarget.randoms.randoms_in_a_brick_from_edges()`.

    Returns
    -------
    :class:`~numpy.ndarray`
        a random catalog with the same columns as returned by
        :func:`~desitwo.randoms.get_quantities_in_a_brick()`
    :class:`str`
        Name of the filename used/passed to read the bright star mask.
    """
    # ADM grab brick information for this directory.
    brickdict = get_brick_info(direc)
    # ADM also grab the (unique) brick names from the dictionary.
    bricknames = np.array(list(brickdict.keys()))

    # ADM a little more information if we're slurming across nodes.
    if os.getenv("SLURMD_NODENAME") is not None:
        log.info(f"Running on Node {os.getenv('SLURMD_NODENAME')}")

    # ADM recover the pixel-level quantities in the DR bricks.
    randoms = select_randoms_bricks(
        brickdict, bricknames, direc, density=density, numproc=numproc,
        dustdir=dustdir, aprad=aprad, seed=seed)

    # ADM retrieve whether the random location was in a bright star mask.
    # ADM there's no point parallelizing this, it's very fast.
    inMx, Mxfn = is_in_gaia_mask(randoms["RA"], randoms["DEC"])
    # ADM add whether the location is in a bright star mask.
    randoms["IN_ARJUN_MASK"] = inMx

    # ADM one last shuffle to randomize across brick boundaries.
    np.random.seed(615+seed)
    np.random.shuffle(randoms)

    return randoms, Mxfn
