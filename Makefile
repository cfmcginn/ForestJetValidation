#std-c++11
CXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs  -Werror -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ROOT=`root-config --cflags --glibs`
INCLUDE=-I $(PWD)
MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_PDFDIR=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirPdfDir forestJetValidation plotAllTH2 plotPFEta

mkdirBin:
	$(MKDIR_BIN)

mkdirPdfDir:
	$(MKDIR_PDFDIR)

forestJetValidation: src/forestJetValidation.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/forestJetValidation.exe src/forestJetValidation.C

plotAllTH2: src/plotAllTH2.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/plotAllTH2.exe src/plotAllTH2.C

plotPFEta: src/plotPFEta.C
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ROOT) -o bin/plotPFEta.exe src/plotPFEta.C

clean:
	rm -f *~
	rm -f include/*~
	rm -f include/#*#
	rm -f src/*~
	rm -f src/#*#
	rm -f bin/*.exe
	rmdir bin || true