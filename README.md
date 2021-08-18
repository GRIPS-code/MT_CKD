# Introduction

**NOTE: These versions are a work in progress! Please see the
[original MT_CKD reposity](https://github.com/AER-RC/MT_CKD) for the official version.**

This repository contains code for the AER continuum model, taken from
[LBLRTM](https://github.com/AER-RC/LBLRTM).  The code has been rewritten as a standalone
library.  MT_CKD is an enhancement ot the original CKD model.  Separate FORTRAN-90 and python
versions are provided.  The FORTRAN-90 sources are located in the `src` directory and
can be built using the instructions below.  The python version is located in the `mt_ckd`
directory and can be installed as described below.

# Building the FORTRAN-90 library.

Building this model requires autotools, a fortran compiler, and the netCDF library.
To build, run:

```bash
$ autoreconf --install
$ ./configure FC=<path to fortran compiler> FCFLAGS=-I<path to netcdf includedir> \
              LDFLAGS=-L<path to netcdf libdir>
$ make
```

# Run the FORTRAN-90 Example

To run a simple test, try:

```bash
$ make check
```

This produces an output netCDF dataset `tests/out.nc` containing the continua.

# Python version
A work-in-progess python version is being attempted.  To install the python package, run the
following in the base of the repository.

```bash
$ pip install --upgrade pip
$ pip install .
```

A simple example is provided in the `tests` directory.  Simply run:

```bash
$ python3 tests/test_mt_ckd.py
```
