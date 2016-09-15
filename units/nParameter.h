/**===========================================================================

    lossyWAV: Added noise WAV bit reduction method by David Robinson;
              Noise shaping coefficients by Sebastian Gesemann;

    Copyright (C) 2007-2016 Nick Currie, Copyleft.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: lossywav <at> hotmail <dot> co <dot> uk

==============================================================================
    Initial translation to C++ from Delphi
    Copyright (C) Tyge Løvset (tycho), Aug. 2012
===========================================================================**/

#ifndef nParameter_h_
#define nParameter_h_

static const char lossyWAVHead1 [] = ", Copyright (C) 2007-2016 Nick Currie. Copyleft.\n";
static const char lossyWAVHead2 [] = "This is free software under the GNU GPLv3+ license; There is NO WARRANTY, to\n"
                                     "the extent permitted by law. <http://www.gnu.org/licenses/> for details.\n";
extern std::ofstream LogOutput;

void nParameter_Init(int32_t argc, char* argv[]);

void nParameter_Cleanup();

std::string WAVFilePrintName();
void open_log_file();
void close_log_file();

#endif // nParameter_h_
