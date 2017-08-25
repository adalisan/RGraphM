## Test environments


## R CMD check results

Solaris 32-bit
Not sure, if tests will run, because of RcppGSL Solaris problem

OSX release
Crashes upon running test (empty matrix problem?)

Windows 7 x64 R 3.3.1
* using R version 3.3.1 (2016-06-21)
* using platform: x86_64-w64-mingw32 (64-bit)

* Windows 10 x64 , R 3.4.0
* For Windows 7 x64
Status:  0 ERROR, 0 WARNINGs, 0 NOTEs

* For Linux 
* Ubuntu 12.04.5
r stable
Status: 0 ERROR, 0 WARNINGs, 0 NOTEs

r oldrelease : 0 ERROR, 0 WARNINGs, 0 NOTEs

r-devel 2017-06-07 r72771: : 0 ERROR, 0 WARNINGs, 0 NOTEs

This package assumes that the GSL library is installed and the environment variable LIB_GSL points to its location. For windows, pre-built libraries are available at [Rtools extras](https://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/)
follow the instructions in this [stackoverflow comment](https://stackoverflow.com/a/23666023/394963)

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

  
* FAILURE SUMMARY

