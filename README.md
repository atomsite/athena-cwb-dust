athena
======
<!-- Jenkins Status Badge in Markdown (with view), unprotected, flat style -->
<!-- In general, need to be on Princeton VPN, logged into Princeton CAS, with ViewStatus access to Jenkins instance to click on unprotected Build Status Badge, but server is configured to whitelist GitHub -->
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

<!--[![Public GitHub  issues](https://img.shields.io/github/issues/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/issues)
[![Public GitHub pull requests](https://img.shields.io/github/issues-pr/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/pulls) -->

Athena++ radiation MHD code

This version has been added upon to simulate dusty colliding wind binary systems.

## Compile Instructions

- A working parallel version of HDF5 is a must!
- Configuring using python programme, my configuration is:

```bash
./configure.py -omp -hdf5 --hdf5_path=$HDF5_LOC --prob=cwb --nscalars=4
```

- Where the variable `$HDF5_LOC` is the directory your HDF5 libraries are located

## Notes

- There appears to be an issue with OpenMPI, at least for me, OpenMP or the Intel MPI libraries are recommended instead!