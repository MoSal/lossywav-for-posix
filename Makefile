HEADERS = version.h \
          units/fftw_interface.h \
          units/nComplex.h \
          units/nCore.h \
          units/nFFT.h \
          units/nFillFFT.h \
          units/nInitialise.h \
          units/nMaths.h \
          units/nOutput.h \
          units/nParameter.h \
          units/nProcess.h \
          units/nRemoveBits.h \
          units/nSGNS.h \
          units/nShiftBlocks.h \
          units/nSpreading.h \
          units/nSupport.h \
          units/nWav.h

OBJS = units/fftw_interface.o \
       units/nCore.o \
       units/nFFT.o \
       units/nFillFFT.o \
       units/nInitialise.o \
       units/nOutput.o \
       units/nParameter.o \
       units/nProcess.o \
       units/nRemoveBits.o \
       units/nSGNS.o \
       units/nShiftBlocks.o \
       units/nSpreading.o \
       units/nWav.o \
       lossyWAV.o

COMMON_CXXFLAGS = -std=c++11 -O2 -pipe
DEFINES = -DHAVE_STD_CHRONO_STEADY_CLOCK_NOW -DHAVE_SETPRIORITY -DHAVE_STAT -DHAVE_CHMOD -DHAVE_NANOSLEEP


all: prep $(OBJS) link

prep:
	$(eval override CXXFLAGS = ${COMMON_CXXFLAGS} ${CXXFLAGS} ${DEFINES})

*.o: ${@:.o=.cpp} $(HEADERS)
	${CXX:-g++} -c ${@:.o=.cpp} -o ${@} ${CXXFLAGS}

link: $(OBJS)
	${CXX} ${OBJS} -o lossywav ${LDFLAGS}

clean:
	-rm -f $(OBJS) lossywav
