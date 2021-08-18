lib_LTLIBRARIES = libmtckd.la
bin_PROGRAMS = continuum
TESTS = tests/run-test.sh


libmtckd_la_SOURCES = $(top_srcdir)/src/utils.f90 \
                      $(top_srcdir)/src/co2_continuum.f90 \
                      $(top_srcdir)/src/h2o_continuum.f90 \
                      $(top_srcdir)/src/n2_continuum.f90 \
                      $(top_srcdir)/src/o2_continuum.f90 \
                      $(top_srcdir)/src/o3_continuum.f90 \
                      $(top_srcdir)/src/continuum.f90
libmtckd_la_LIBADD = -lnetcdff
continuum_SOURCES = $(top_srcdir)/src/main.f90
continuum_LDADD = libmtckd.la

$(top_srcdir)/src/co2_continuum.lo: $(top_srcdir)/src/co2_continuum.f90 \
                                    $(top_srcdir)/src/utils.lo

$(top_srcdir)/src/h2o_continuum.lo: $(top_srcdir)/src/h2o_continuum.f90 \
                                    $(top_srcdir)/src/utils.lo

$(top_srcdir)/src/n2_continuum.lo: $(top_srcdir)/src/n2_continuum.f90 \
                                   $(top_srcdir)/src/utils.lo

$(top_srcdir)/src/o2_continuum.lo: $(top_srcdir)/src/o2_continuum.f90 \
                                   $(top_srcdir)/src/utils.lo

$(top_srcdir)/src/o3_continuum.lo: $(top_srcdir)/src/o3_continuum.f90 \
                                   $(top_srcdir)/src/utils.lo

$(top_srcdir)/src/continuum.lo: $(top_srcdir)/src/continuum.f90 \
                                $(top_srcdir)/src/utils.lo \
                                $(top_srcdir)/src/co2_continuum.lo \
                                $(top_srcdir)/src/h2o_continuum.lo \
                                $(top_srcdir)/src/n2_continuum.lo \
                                $(top_srcdir)/src/o2_continuum.lo \
                                $(top_srcdir)/src/o3_continuum.lo

$(top_srcdir)/src/main.$(OBJEXT): $(top_srcdir)/src/main.f90 \
                                  $(top_srcdir)/src/continuum.lo \
                                  $(top_srcdir)/src/utils.lo

CLEANFILES = *.mod
