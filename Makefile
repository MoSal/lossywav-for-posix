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

CXX := g++ -std=c++11

COMMON_CXXFLAGS = -pipe -Wall -Wextra -O2
COMMON_LDFLAGS = -Wl,-O1,--sort-common,--as-needed
LTO_FLAGS := -flto

OPTIMIZED_CXXFLAGS := -march=native -Ofast

all: prep $(OBJS) link
optimized: prep-optimized $(OBJS) link
fftw: prep-fftw $(OBJS) link
fftw-optimized: prep-fftw-optimized $(OBJS) link

prep:
	$(eval override CXXFLAGS = ${COMMON_CXXFLAGS} ${CXXFLAGS})
	$(eval override LDFLAGS = ${COMMON_LDFLAGS} ${LDFLAGS})
	$(eval override LIBS = ${COMMON_LIBS} ${LIBS})

prep-optimized: prep
	$(eval override CXXFLAGS += ${OPTIMIZED_CXXFLAGS} ${LTO_FLAGS})
	$(eval override LDFLAGS += ${LTO_FLAGS})

prep-fftw: prep
	$(eval override CXXFLAGS += -DENABLE_FFTW)
	$(eval override LIBS += -lfftw3)

prep-fftw-optimized: prep-optimized prep-fftw

*.o: ${@:.o=.cpp} $(HEADERS)
	${CXX} -c ${@:.o=.cpp} -o ${@} ${CXXFLAGS}

link: $(OBJS)
	${CXX} ${OBJS} -o lossywav ${LDFLAGS} ${LIBS}

clean:
	-rm -f $(OBJS) lossywav
