## Test environments
* CentOS 6 , R 3.2.3
* Windows 10 , R 3.2.3
## R CMD check results

For Windows
Status: 1 ERROR, 3 WARNINGs, 1 NOTE

ERROR is due to inconsolata.sty not existing in MikTex distribution in Windows. This is a known issue and there are work arounds

For Linux 
Status: 0 ERROR, ? WARNINGs, ? NOTE


This package assumes that the GSL library is installed and the environment variable LIB_GSL points to its location. For windows, pre-built libraries are available at [Rtools extras](https://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/)
follow the instructions in this [stackoverflow comment](https://stackoverflow.com/a/23666023/394963)

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

  
* FAILURE SUMMARY

