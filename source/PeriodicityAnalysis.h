//
//  PeriodicityAnalysis.h
//  Multipitch
//
//  Created by Nicholas Collins on 08/07/2011.
//  Copyright 2011 Nick Collins. All rights reserved.
//
//this code is under GNU GPL3 license, see file COPYING enclosed

#include "MultipitchHeader.h"

class PeriodicityAnalysis {

public:
    //polyphony estimation iteration data
    float * Ur_; 
    float * Ud_; 
    int maxvoices_; 
    int numvoicesdetected_; 
    float prevmixturescore_; 
    float mixturescore_; 
    float * voiceperiods_; 
    float * voicesaliences_; 
    
    //most salient period estimation 
    int maxq_, q_, qbest_; 
    float *smax_;
    float *torlow_, *torup_;  
    float torprec_, tormin_, tormax_; 
    
    //for results
    float winningtor_; 
    float bestsalience_; 

    float K_; 
    
    float cancellationweight_; 
    float polyphonyestimategamma_;
    
    
    PeriodicityAnalysis(int maxvoices=4); 
    ~PeriodicityAnalysis(); 
    
    //for one frame
    void compute(float * Uk); 
    void minSearch(float * Ur); 
    float smax(int q, float * Ur); 
    
};