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
  Pi0EventDisplayDump.h

lib_LTLIBRARIES = \
  libPi0EventDisplayDump.la

libPi0EventDisplayDump_la_SOURCES = \
  Pi0EventDisplayDump.cc

libPi0EventDisplayDump_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lcalo_io \
  -lg4dst \
  -lphg4hit

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libPi0EventDisplayDump.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
