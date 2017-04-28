## Test environments
* CentOS 6 , R 3.2.3
* Windows 10 , R 3.2.3
## R CMD check results

For Windows
Status: 1 ERROR, 3 WARNINGs, 1 NOTE

ERROR is due to inconsolata.sty not existing in MikTex distribution in Windows. This is a known issue and there are work arounds

For Linux 
Status: 0 ERROR, 0 WARNINGs, 4 NOTEs
Undocumented code objects:
  ‘a’ ‘b’
Undocumented data sets:
  ‘a’ ‘b’
Package has both ‘src/Makevars.in’ and ‘src/Makevars’.
Installation with --no-configure' is unlikely to work.  If you intended
‘src/Makevars’ to be used on Windows, rename it to ‘src/Makevars.win’
otherwise remove it.  If ‘configure’ created ‘src/Makevars’, you need a
‘cleanup’ script.
* checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS) ... OK
* checking compiled code ... NOTE
File ‘/mnt/k/Sancar/Documents/projects/rgraphm_linux.Rcheck/RGraphM/libs/RGraphM.so’:
  Found ‘puts’, possibly from ‘printf’ (C), ‘puts’ (C)
    Objects: ‘./graphm/./experiment.o’, ‘./graphm/./graph.o’

Compiled code should not call entry points which might terminate R nor
write to stdout/stderr instead of to the console.


This package assumes that the GSL library is installed and the environment variable LIB_GSL points to its location. For windows, pre-built libraries are available at [Rtools extras](https://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/)
follow the instructions in this [stackoverflow comment](https://stackoverflow.com/a/23666023/394963)

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

  
* FAILURE SUMMARY

