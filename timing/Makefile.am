
## Makefile.am -- Process this file with automake to produce Makefile.in

## compile library

noinst_PROGRAMS = \
    ddedge \
    fids \
    mark \
    opids \
    reduce \
    unpacked

AM_CXXFLAGS = -Wall

ddedge_SOURCES =    timer.h \
                    park_random.h \
                    reporting.h reporting.cc \
                    ddedge.cc
ddedge_LDADD =      ../src/libmeddly.la

#####################################################################

fids_SOURCES =      timer.h \
                    park_random.h \
                    reporting.h reporting.cc \
                    fids.cc
fids_LDADD =        ../src/libmeddly.la

#####################################################################

mark_SOURCES =      timer.h \
                    reporting.h reporting.cc \
                    simple_model.h simple_model.cc \
                    mark.cc
mark_LDADD =        ../src/libmeddly.la

#####################################################################

reduce_SOURCES =    timer.h \
                    park_random.h \
                    reporting.h reporting.cc \
                    reduce.cc
reduce_LDADD =      ../src/libmeddly.la

#####################################################################

opids_SOURCES =     timer.h \
                    park_random.h \
                    reporting.h reporting.cc \
                    opids.cc
opids_LDADD =       ../src/libmeddly.la

#####################################################################

unpacked_SOURCES =  timer.h \
                    park_random.h \
                    reporting.h reporting.cc \
                    unpacked.cc
unpacked_LDADD =    ../src/libmeddly.la

