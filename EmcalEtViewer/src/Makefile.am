AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  DisplayLegoPlot.h 

lib_LTLIBRARIES = \
  libDisplayLegoPlot.la

libDisplayLegoPlot_la_SOURCES = \
  DisplayLegoPlot.cc

libDisplayLegoPlot_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lffaobjects \
  -lffarawobjects \
  -lcalo_io \
  -ltrack_io \
  -lintt

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libDisplayLegoPlot.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
