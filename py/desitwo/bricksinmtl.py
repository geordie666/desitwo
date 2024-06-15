"""
desitwo.bricksinmtl
===================

Determine the set of bricks covered by a set of dark/bright MTL ledgers.

"""
import os
import numpy as np
import matplotlib.pyplot as plt
from time import time
from glob import glob

from desitarget.mtl import get_mtl_dir
from desitarget import io

from desiutil import brick
from desiutil.log import get_logger
from desiutil.plots import plot_sky_binned

# ADM set up the Legacy Surveys bricks object.
bricks = brick.Bricks(bricksize=0.25)

# ADM set up the DESI default logger.
log = get_logger()


def bricks_in_mtl(mtldir=None, plotfile=None, onlyobserved=False):
    """Determine set of bricks that touches bright/dark MTL ledgers.

    Parameters
    ----------
    mtldir : :class:`str`, optional, defaults to ``None``
        Full path to the directory that hosts the MTL ledgers and the MTL
        tile file. If ``None``, then look up the MTL directory from the
        $MTL_DIR environment variable.
    plotfile : :class:`str`, optional, defaults to ``None``
        Full path to file to which to write a plot of the locations of
        the bricks. If ``None`` then do not produce a plot.
    onlyobserved : :class:`bool`, optional, defaults to ``False``
        If ``True`` only return bricks touched by targets that have
        actually been observed.

    Returns
    -------
    :class:`~numpy.array`
        The standard Table of brick information, for bricks that touch
        the ledgers.
    """
    # ADM start the clock.
    t0 = time()

    # ADM grab MTL directory, returning environment variable for None.
    mtldir = get_mtl_dir(mtldir)

    # ADM to hold the list of brick IDs.
    bids = []
    # ADM loop over both the bright and dark ledgers.
    for obscon in "bright", "dark":
        # ADM retrieve all of the ledgers for bright/dark.
        mdir = io.find_target_files(mtldir, flavor="mtl", obscon=obscon)
        fns = sorted(glob(os.path.join(mdir, "mtl*ecsv")))
        # ADM loop over all the ledgers...
        for i, fn in enumerate(fns):
            # ADM ...read in the ledgers...
            mtl = io.read_mtl_ledger(fn)
            # ADM if requested, limit to just observed targets.
            if onlyobserved:
                ii = mtl["NUMOBS"] > 0
                mtl = mtl[ii]
            # ADM ...and determine the bricks covered by the ledgers.
            bids += list(bricks.brickid(mtl["RA"], mtl["DEC"]))
            # ADM restrict to only unique bricks.
            bids = list(set(bids))
            if i % 500 == 499:
                msg = f"{i+1}/{len(fns)} {obscon} bricks...t={time()-t0:.1f}s"
                log.info(msg)

    # ADM retrieve the full brick information for the relevant bricks.
    bids = np.array(sorted(bids))
    bricktable = bricks.to_table()

    # ADM remember the -1 as the table is zero-indexed but
    # ADM starts at BRICKID==1
    bricktable = bricktable[bids-1]

    # ADM produce a plot, if requested.
    if plotfile is not None:
        plot_sky_binned(bricktable["RA"], bricktable["DEC"])
        plt.savefig(plotfile)
        msg = f"Plot written to {plotfile}...t={time()-t0:.1f}s"
        log.info(msg)

    return bricktable
