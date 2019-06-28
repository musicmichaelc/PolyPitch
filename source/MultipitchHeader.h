//
//  MultipitchHeader.h
//  Multipitch
//
//  Created by Nicholas Collins on 29/05/2011.

//this code is under GNU GPL3 license, see file COPYING enclosed

#pragma once

#include <iostream>
#include <math.h>
#include "samplerate.h"
#include "SC_PlugIn.h"

//sampling rate assumed at 44100
//was const float, now global and can be changed 
extern float g_samplingrate; 
extern float g_sampleperiod; 

const int g_ppwindowsize = 2048; 
const int g_ppdecimatedwindowsize = 256; 
const int g_Uksize = 256; 

extern InterfaceTable *ft; 


//build with accelerate framework directly, bypassing SC fft (new version of which is SC 3.5+ only)
#define OSXACCELERATENOTSCFFT