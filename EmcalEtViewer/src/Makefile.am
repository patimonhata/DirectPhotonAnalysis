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
  CemcTowerDumper.h \
  CemcClusterDumper.h 

lib_LTLIBRARIES = \
  libCemcTowerDumper.la \
  libCemcClusterDumper.la

libCemcTowerDumper_la_SOURCES = \
  CemcTowerDumper.cc

libCemcTowerDumper_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lffaobjects \
  -lffarawobjects \
  -lcalo_io \
  -ltrack_io \
  -lintt

libCemcClusterDumper_la_SOURCES = \
  CemcClusterDumper.cc

libCemcClusterDumper_la_LIBADD = \
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
testexternals_LDADD   = libCemcTowerDumper.la libCemcClusterDumper.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
