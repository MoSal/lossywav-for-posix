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

COMMON_CXXFLAGS = -std=c++11 -O2 -pipe -Wall -Wextra
DEFINES = -DHAVE_STD_CHRONO_STEADY_CLOCK_NOW -DHAVE_SETPRIORITY -DHAVE_STAT -DHAVE_CHMOD
COMMON_LDFLAGS = -Wl,-O1 -Wl,--sort-common -Wl,--as-needed
LTO_FLAGS := -flto

OPTIMIZED_CXXFLAGS := -march=native -Ofast

all: prep $(OBJS) link
optimized: prep-optimized $(OBJS) link
fftw: prep-fftw $(OBJS) link
fftw-optimized: prep-fftw-optimized $(OBJS) link

prep:
	$(eval override CXXFLAGS = ${COMMON_CXXFLAGS} ${CXXFLAGS} ${DEFINES})
	$(eval override LDFLAGS = ${COMMON_LDFLAGS} ${LDFLAGS})
	$(eval override LIBS = ${COMMON_LIBS} ${LIBS})

prep-optimized: prep
	$(eval override CXXFLAGS += ${OPTIMIZED_CXXFLAGS} ${LTO_FLAGS})
	$(eval override LDFLAGS += ${LTO_FLAGS})

prep-fftw: prep
	$(eval override CXXFLAGS += -DHAVE_FFTW3)
	$(eval override LIBS += -lfftw3)

prep-fftw-optimized: prep-optimized prep-fftw

*.o: ${@:.o=.cpp} $(HEADERS)
	${CXX:-g++} -c ${@:.o=.cpp} -o ${@} ${CXXFLAGS}

link: $(OBJS)
	${CXX} ${OBJS} -o lossywav ${LDFLAGS} ${LIBS}

clean:
	-rm -f $(OBJS) lossywav
