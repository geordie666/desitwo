.. title: Focal Plane
.. slug: focalplane
.. tags: 
.. has_math: yes

.. |sigma|    unicode:: U+003C3 .. GREEK SMALL LETTER SIGMA
.. |sup2|     unicode:: U+000B2 .. SUPERSCRIPT TWO
.. |alpha|      unicode:: U+003B1 .. GREEK SMALL LETTER ALPHA
.. |chi|      unicode:: U+003C7 .. GREEK SMALL LETTER CHI
.. |delta|    unicode:: U+003B4 .. GREEK SMALL LETTER DELTA
.. |deg|    unicode:: U+000B0 .. DEGREE SIGN
.. |times|  unicode:: U+000D7 .. MULTIPLICATION SIGN
.. |plusmn| unicode:: U+000B1 .. PLUS-MINUS SIGN
.. |Prime|    unicode:: U+02033 .. DOUBLE PRIME
.. |geq|    unicode:: U+02265 .. GREATER THAN OR EQUAL TO


.. class:: pull-right well

.. contents::

The focal plane is composed of a single mechanical structure that accepts the
fiber positioners, the guide/focus cameras, fiducial reference sources, and sky
monitor fibers. The fiber positioners are modularized in 204 “rafts”. Each raft
is an equilateral triangle measuring 74 mm on a side and 604 mm long. The raft
is a self-contained unit with 63 positioners and control electronics. One fiber
cable and three Power-over-Ethernet (PoE) cables attach at the rear. The present
design weighs only 1.2 kg per raft, making them easily removable and serviceable.

.. container:: col-md-5 col-left

  Fiber positioning robots are mounted at a 6.2 mm center-to-center pitch. The raft
  module can accommodate any of several fiber positioning technologies, so long as they can
  be mounted at 6.2 mm pitch. Our *Reference Design* for the raft incorporates 21 “Trillium”
  robot units. Each Trillium has 3 fiber actuators, driven by 6 brushless DC gear motors of
  diameter 4 mm, these are the identical motors that drive the DESI robots. The change from the
  DESI design that allows for a reduction in pitch was the introduction of an additional set of
  transfer gears which couple the eccentric axis to the central rotation axis. This allows the
  motors to be axially nested, such that one motor is mounted partially behind the other.

.. container:: col-md-7 col-right

   .. image:: ../../files/instrument/focalplane.jpg
      :height: 375
      :width: 800



