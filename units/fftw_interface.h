/**------------------------------------------------------------
  Delphi interface to the FFTW library -- FFTW Version 3.0.1.
  Note that this interface is incomplete. Additional function
  interface entries may be added in an anologous manner, see
  FFTW for more details.

  Last modified 22/DEC/03
  Written by Mark G. Beckett (g.beckett@epcc.ed.ac.uk
 --------------------------------------------------------------
  Modified 07/JUN/2009
  by Nick Currie (lossyWAV at yahoo dot com)
  Now dynamically links rather than statically. Required to
  stop program crashing if libFFTW3-3.DLL is not available.
 --------------------------------------------------------------
  Modified 14/APR/2010
  by Nick Currie (lossyWAV at hotmail <dot> co <dot> uk
 --------------------------------------------------------------
    Initial translation to C++ from Delphi
    Copyright (C) Tyge Løvset (tycho), Aug. 2012
 -----------------------------------------------------------**/

#ifndef fftw_interface_h_
#define fftw_interface_h_

#include "nCore.h"

struct FFTW_Rec
{
    void* (* Plan_DFT_r2c_1d)(int32_t n, double* inData, double* outData, uint32_t flags);
    void* (* Plan_DFT_1d)(int32_t n, double* inData, double* outData, int32_t Sign, uint32_t flags);
    void (* Execute_R2C_New_Array)(void* plan, double* inData, double* outData);
    void (* Execute_C2C_New_Array)(void* plan, double* inData, double* outData);
    void (* Destroy_Plan)(void* plan);
    void* Plans [MAX_FFT_BIT_LENGTH + 1];
    void* Plans_Inv [MAX_FFT_BIT_LENGTH + 1];
};

enum
{
    FFTW_FORWARD = -1,
    FFTW_BACKWARD = 1,
    FFTW_MEASURE = 0,
    FFTW_DESTROY_INPUT      = 1,
    FFTW_UNALIGNED          = 1 << 1,
    FFTW_CONSERVE_MEMORY    = 1 << 2,
    FFTW_EXHAUSTIVE         = 1 << 3,
    FFTW_PRESERVE_INPUT     = 1 << 4,
    FFTW_PATIENT            = 1 << 5,
    FFTW_ESTIMATE           = 1 << 6
};

extern FFTW_Rec FFTW;

bool FFTW_Initialised();

bool FFTW_Initialise();

void FFTW_Cleanup();

#endif // fftw_interface_h_
