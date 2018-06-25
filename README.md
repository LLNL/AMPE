AMPE v2.0
=========

Adaptive Mesh Phase-field Evolution

Authors
-------

 * Milo R. Dorr (dorr1@llnl.gov)
 * Jean-Luc Fattebert (fattebertj@ornl.gov)
 * Mike E. Wickett (wickett1@llnl.gov)

Dependencies
------------

* [Hypre] (https://github.com/LLNL/hypre)

* [SAMRAI] (https://github.com/LLNL/SAMRAI)

* [CPODES] (https://simtk.org/projects/cpodes)

* [Sundials] (https://github.com/LLNL/sundials)

* [HDF5] (https://support.hdfgroup.org/HDF5)

* [NetCDF] (https://www.unidata.ucar.edu/software/netcdf)

* [BOOST] (http://www.boost.org)

Since CPODES is currently not distributed with Sundials, and SAMRAI
does not supports the lastest SUNDIALS release, modifications to 
these libraries had to be made. The modified libraries are distributed
with this code under the 'base' directory and need to be built before
building the main code.

References
----------

M. R. Dorr, J.-L. Fattebert, M. E. Wickett, J. F. Belak, P. E. A. Turchi,
"A Numerical Algorithm for the Solution of a Phase-Field Model of
Polycrystalline materials",
J. Comp. Phys. 229 (3), p. 626-641 (2010)

J.-L. Fattebert, M. E. Wickett, P. E. A. Turchi, 
"Phase-field modeling of coring during solidification of Au-Ni alloy using
quaternions and CALPHAD input",
ACTA MATERIALIA, 62, (2014), 89-104

Release
-------

Copyright (c) 2018, Lawrence Livermore National Security, LLC
and UT-Battelle, LLC.

Produced at the Lawrence Livermore National Laboratory and
the Oak Ridge National Laboratory.

All rights reserved.

Unlimited Open Source - BSD Distribution

For release details and restrictions, please read the LICENSE file.
It is also linked here:
- [LICENSE](./LICENSE)

`LLNL-CODE-747500`  `OCEC-18-028`

