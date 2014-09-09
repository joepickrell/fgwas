
bin_PROGRAMS = fgwas test
DISTCHECK_CONFIGURE_FLAGS=LDFLAGS=-L/opt/local/lib CXXFLAGS=-I/opt/local/include
fgwas_SOURCES = CmdLine.cpp fgwas.cpp gzstream.cpp SNP.cpp SNPs.cpp fgwas_params.cpp
test_SOURCES = test.cpp CmdLine.cpp gzstream.cpp SNP.cpp SNPs.cpp fgwas_params.cpp

noinst_HEADERS = gzstream.h CmdLine.h SNP.h SNPs.h fgwas_params.h
