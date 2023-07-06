#  Optional configure Makefile overrides for bcftools.
#
#    Copyright (C) 2015,2017 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# This is config.mk.  Generated from config.mk.in by configure.
#
# If you use configure, this file overrides variables and augments rules
# in the Makefile to reflect your configuration choices.  If you don't run
# configure, the main Makefile contains suitable conservative defaults.

prefix       = /usr/local
exec_prefix  = ${prefix}
bindir       = ${exec_prefix}/bin
libexecdir   = ${exec_prefix}/libexec
datarootdir  = ${prefix}/share
mandir       = ${datarootdir}/man

CC       = gcc
CPPFLAGS = 
CFLAGS   =  -Wall -g -O2
LDFLAGS  = 
LIBS     = -ldl 

DYNAMIC_FLAGS = -rdynamic

USE_GPL = 
GSL_LIBS = 
PERL_CFLAGS = 
PERL_LIBS = 

PLATFORM   = default
PLUGINS_ENABLED = yes
plugindir  = $(libexecdir)/bcftools
pluginpath = $(libexecdir)/bcftools
PLUGIN_EXT = .so

HTSDIR = htslib-1.11
include $(HTSDIR)/htslib.mk
include $(HTSDIR)/htslib_static.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSLIB_LIB = $(HTSLIB) $(HTSLIB_static_LIBS)
HTSLIB_DLL = $(HTSDIR)/
HTSLIB_LDFLAGS = $(HTSLIB_static_LDFLAGS)
W32_PLUGIN_LIBS = libbcftools.a $(HTSLIB_DLL) $(ALL_LIBS)
BGZIP = $(HTSDIR)/bgzip
TABIX = $(HTSDIR)/tabix
HTSLIB_CPPFLAGS = -Ihtslib-1.11
#HTSLIB_LDFLAGS = -Lhtslib-1.11
#HTSLIB_LIB = -lhts
#W32_PLUGIN_LIBS = libbcftools.a $(HTSLIB_LDFLAGS) $(HTSLIB_LIB) $(ALL_LIBS)
