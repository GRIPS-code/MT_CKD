lib_LTLIBRARIES = libmtckd.la libmtckd2.la
bin_PROGRAMS = orig continuum
TESTS = tests/run-test.sh
libmtckd_la_SOURCES = $(top_srcdir)/src/contnm.f90 \
                      $(top_srcdir)/src/lblparams.f90 \
                      $(top_srcdir)/src/phys_consts.f90 \
                      $(top_srcdir)/src/planet_earth.f90
$(top_srcdir)/src/contnm.lo: $(top_srcdir)/src/contnm.f90 \
                             $(top_srcdir)/src/lblparams.lo \
                             $(top_srcdir)/src/phys_consts.lo
$(top_srcdir)/src/planet_earth.lo: $(top_srcdir)/src/planet_earth.f90 \
                                   $(top_srcdir)/src/phys_consts.lo
orig_SOURCES = $(top_srcdir)/src/cntnm_progr.f90
orig_LDADD = libmtckd.la


libmtckd2_la_SOURCES = $(top_srcdir)/src/utils.f90 \
                       $(top_srcdir)/src/co2_continuum.f90 \
                       $(top_srcdir)/src/h2o_continuum.f90 \
                       $(top_srcdir)/src/n2_continuum.f90 \
                       $(top_srcdir)/src/o2_continuum.f90 \
                       $(top_srcdir)/src/o3_continuum.f90 \
                       $(top_srcdir)/src/continuum.f90
continuum_SOURCES = $(top_srcdir)/src/main.f90
continuum_LDADD = libmtckd2.la -lnetcdff


CLEANFILES = *.mod *__genmod.f90
