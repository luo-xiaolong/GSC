/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if BCFtools should enable plugins. */
#define ENABLE_BCF_PLUGINS 1

/* Define if BCFtools should enable for support PERL scripts in -i/-e
   filtering expressions. */
/* #undef ENABLE_PERL_FILTERS */

/* Define to 1 if you have the `z' library (-lz). */
#define HAVE_LIBZ 1

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "samtools-help@lists.sourceforge.net"

/* Define to the full name of this package. */
#define PACKAGE_NAME "BCFtools"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "BCFtools 1.11"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "bcftools"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://www.htslib.org/"

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.11"

/* Platform-dependent plugin filename extension. */
#define PLUGIN_EXT ".so"

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */
