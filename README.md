SkyNet
======

SkyNet is a general-purpose nuclear reaction network for nuclear astrophysics
applications.

SkyNet is available on github: https://bitbucket.org/jlippuner/skynet

The original author and primary maintainer of SkyNet is Jonas Lippuner
(jonas@lippuner.ca). Luke Roberts (robertsl@nscl.msu.edu) has made
substantial contributions.

A detailed methods paper describing SkyNet and the physics it implements has
been published:

  Lippuner, J. and L. F. Roberts (2017), "SkyNet: A modular nuclear reaction
  network library", submitted to ApJS, https://arxiv.org/abs/1706.06198

The above paper should be cited in any academic work that uses SkyNet or is
based on SkyNet.



DOCUMENTATION
=============

It should be possible to compile and run SkyNet with any major Linux
distribution and Mac OS. Microsoft Windows is not supported. See the INSTALL
file for details about how to compile and run SkyNet.

Examples are provided in the examples directory.

For instructions on how to make movies with SkyNet, see README_make_movie.

There are a number of README files for specific machines where SkyNet has been
tested.



LICENSE
=======

SkyNet is distributed under the Revised BSD 3-Clause License. See the LICENSE
file for details.



ACKNOWLEDGEMENTS
================

SkyNet can read JINA REACLIB database files that contain reaction rate fits for
various nuclear reactions. A copy of a REACLIB snapshot is distributed with
SkyNet in the data directory. REACLIB is available at
https://groups.nscl.msu.edu/jina/reaclib/db/index.php

SkyNet uses TinyXML-2 developed by Lee Thomason, which is available at
https://github.com/leethomason/tinyxml2

SkyNet contains the Helmholtz equation of state developed by Frank Timmes, which
is available at http://cococubed.asu.edu/code_pages/eos.shtml (specifically
http://cococubed.asu.edu/codes/eos/helmholtz.tbz). The Helmhotz EOS has been
described in:

  Timmes, F. X. and F. D. Swesty (2000), "The Accuracy, Consistency, and Speed
  of an Electron-Positron Equation of State Based on Table Interpolation of
  the Helmholtz Free Energy", Astrophysical Journal, Supplement 126, p. 501,
  doi: 10.1086/313304

The X-ray burst trajectory used in the tests was graciously provided by Hendrik
Schatz. For details, see

  Schatz, H., A. Aprahamian, V. Barnard, L. Bildsten, A. Cumming, M. Ouellette,
  T. Rauscher, F.-K. Thielemann, and M. Wiescher (2001), "End Point of the rp
  Process on Accreting Neutron Stars", Physical Review Letters 86, p. 3471.
  doi: 10.1103/PhysRevLett.86.3471, https://arxiv.org/abs/astro-ph/0102418

The original development of SkyNet has been funded by the United States National
Science Foundation (NSF) under the Theoretical and Computational Astrophysics
Networks (TCAN) program (award number AST-1333520).

This project has also been supported in part by the Sherman Fairchild Foundation
and the JINA Center for the Evolution of the Elements (NSF grant PHY-1430152).
