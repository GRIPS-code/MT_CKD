lib_LTLIBRARIES = libmtckd.la
bin_PROGRAMS = continuum
TESTS = tests/run-test.sh

libmtckd_la_SOURCES = $(top_srcdir)/src/contnm.f90 \
                      $(top_srcdir)/src/lblparams.f90 \
                      $(top_srcdir)/src/phys_consts.f90 \
                      $(top_srcdir)/src/planet_earth.f90

continuum_SOURCES = $(top_srcdir)/src/cntnm_progr.f90
continuum_LDADD = libmtckd.la

$(top_srcdir)/src/contnm.lo: $(top_srcdir)/src/contnm.f90 \
                             $(top_srcdir)/src/lblparams.lo \
                             $(top_srcdir)/src/phys_consts.lo

$(top_srcdir)/src/planet_earth.lo: $(top_srcdir)/src/planet_earth.f90 \
                                   $(top_srcdir)/src/phys_consts.lo

CLEANFILES = *.mod *__genmod.f90
