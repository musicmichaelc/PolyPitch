//
//  PeriodicityAnalysis.h
//  Multipitch
//
//  Created by Nicholas Collins on 19/07/2011.
//  Copyright 2011 Nick Collins. All rights reserved.
//
//this code is under GNU GPL3 license, see file COPYING enclosed

#include "GammatoneFilterBank.h"
#include "PeriodicityAnalysis.h"



class PolyPitchUGen {

public:
    //polyphony estimation iteration data
    
    GammatoneFilterBank gamma_; // = GammatoneFilterBank(70, 65.f, 5200.f);
    PeriodicityAnalysis polypitch_;

    float inputbuffer_[g_ppwindowsize]; 
    float holdbuffer_[g_ppwindowsize]; 
    
    int pos_; 
    
    int amortisationstage_; 
    
    PolyPitchUGen(float samplingRate,int maxvoices, World * world); 
    ~PolyPitchUGen(); 
    
    //numSamples should be 64! 
    void compute(float * input, int numSamples); 
        
};