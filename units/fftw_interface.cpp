// --------------------------------------------------------------
//  Delphi interface to the FFTW library -- FFTW Version 3.0.1.
//  Note that this interface is incomplete. Additional function
//  interface entries may be added in an anologous manner, see
//  FFTW for more details.
//
//  Last modified 22/DEC/03
//  Written by Mark G. Beckett (g.beckett@epcc.ed.ac.uk
// --------------------------------------------------------------
//  Modified 07/JUN/2009
//  by Nick Currie (lossyWAV at yahoo dot com)
//  Now dynamically links rather than statically. Required to
//  stop program crashing if libFFTW3-3.DLL is not available.
// --------------------------------------------------------------
//  Modified 14/APR/2010
//  by Nick Currie (lossyWAV at hotmail <dot> co <dot> uk
// --------------------------------------------------------------
//    Initial translation to C++ from Delphi
//    by Tyge Løvset (tycho), Aug. 2012
// --------------------------------------------------------------

#include "fftw_interface.h"

static HMODULE FFTW_DLL_Handle = nullptr;
static bool FFTW_DLL_Loaded = false;

FFTW_Rec FFTW;

bool Check_Initialised(FFTW_Rec& FFTW_Record)
{
    return ((&FFTW_Record.Plan_DFT_r2c_1d != nullptr)  && (&FFTW_Record.Plan_DFT_1d != nullptr) && (&FFTW_Record.Execute_R2C_New_Array != nullptr) && (&FFTW_Record.Destroy_Plan != nullptr));
}


bool FFTW_Initialised()
{
    return FFTW_DLL_Loaded && (FFTW_DLL_Handle != nullptr);
}


void Initialise_FFTW_DLL(const char* this_DLL)
{
    FFTW_DLL_Handle = LoadLibrary(this_DLL);

    if (FFTW_DLL_Handle != nullptr)
    {
        FFTW.Plan_DFT_r2c_1d = (void*(*)(int,double*,double*,uint32_t)) GetProcAddress(FFTW_DLL_Handle, "fftw_plan_dft_r2c_1d");
        FFTW.Plan_DFT_1d = (void*(*)(int,double*,double*,int,uint32_t)) GetProcAddress(FFTW_DLL_Handle, "fftw_plan_dft_1d");

        FFTW.Execute_R2C_New_Array = (void(*)(void*, double*, double*)) GetProcAddress(FFTW_DLL_Handle, "fftw_execute_dft_r2c");
        FFTW.Execute_C2C_New_Array = (void(*)(void*, double*, double*)) GetProcAddress(FFTW_DLL_Handle, "fftw_execute_dft");

        FFTW.Destroy_Plan = (void(*)(void*)) GetProcAddress(FFTW_DLL_Handle, "fftw_destroy_plan");
    }

    FFTW_DLL_Loaded = Check_Initialised(FFTW) && (FFTW_DLL_Handle != nullptr);
}


bool FFTW_Initialise()
{
    if (FFTW_DLL_Handle==nullptr)
    {
        Initialise_FFTW_DLL((char*) "libfftw3-3.dll\0");
    }

    // Kludge to allow 64-bit DLL to be loaded (if renamed file exists on path) for 64-bit build.
    if (FFTW_DLL_Handle==nullptr)
    {
        Initialise_FFTW_DLL((char*) "libfftw3-3_64.dll\0");
    }

    return FFTW_Initialised();
}


bool CHECK_FFTW_DLL_Loaded()
{
    return FFTW_DLL_Loaded;
}


void FFTW_Cleanup()
{
    if (FFTW_DLL_Handle != nullptr)
    {
        for (int32_t fc_i = 1; fc_i != MAX_FFT_BIT_LENGTH; ++fc_i)
        {
            if (FFTW.Plans[fc_i] != nullptr)
            {
                FFTW.Destroy_Plan(FFTW.Plans[fc_i]);
            }
            if (FFTW.Plans_Inv[fc_i] != nullptr)
            {
                FFTW.Destroy_Plan(FFTW.Plans_Inv[fc_i]);
            }
        }

        FFTW.Plan_DFT_r2c_1d = nullptr;
        FFTW.Plan_DFT_1d = nullptr;
//        FFTW.Execute = nullptr;
        FFTW.Execute_C2C_New_Array = nullptr;
        FFTW.Execute_R2C_New_Array = nullptr;
        FFTW.Destroy_Plan = nullptr;

        FreeLibrary(FFTW_DLL_Handle);
        FFTW_DLL_Handle = nullptr;
        FFTW_DLL_Loaded = false;
    }
}
