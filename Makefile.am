# -*- Automake -*-
#
# Makefile for the AtomPAW package
#
# Copyright (C) 2010 Yann Pouillon
#
# This file is part of the AtomPAW software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Top makefile
#

# Autotools parameters
ACLOCAL_AMFLAGS = -I config/m4

# Subdirectories containing makefiles
SUBDIRS = src

# Files to distribute in the source package
EXTRA_DIST = \
	README \
	config/m4/ax_f90_module_extension.m4 \
	config/m4/fortran.m4 \
	config/m4/information.m4 \
	config/m4/linalg.m4 \
	config/m4/optimizations.m4 \
	config/m4/workarounds.m4
 
                    # ------------------------------------ #

#
# Documentation: examples
#

# Files to distribute in the source package
atp_xmpl_srcs = \
	example/Cu/Cu.atomicdata \
	example/Cu/Cu.atompaw.input \
	example/C/C.atomicdata \
	example/C/C.atompaw.input \
	example/Li/Li.atomicdata.gz \
        example/Li/Li.GGA-PBE-paw.abinit.gz \
        example/Li/Li.GGA-PBE-paw.UPF.gz \
	example/Li/Li.atompaw.input \
	example/F/F.atomicdata.gz \
        example/F/F.GGA-PBE-paw.abinit.gz \
        example/F/F.GGA-PBE-paw.UPF.gz \
	example/F/F.atompaw.input
atp_doc_srcs = \
	doc/atompaw-usersguide.pdf

# Files to install (nobase = preserving the directory structure)
nobase_doc_DATA = $(atp_xmpl_srcs) $(atp_doc_srcs)

# Files to distribute in the source package
EXTRA_DIST += $(atp_xmpl_srcs)
EXTRA_DIST += $(atp_doc_srcs)

