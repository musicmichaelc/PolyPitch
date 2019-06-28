//
//  PeriodicityAnalysis.cpp
//  Multipitch
//
//  Created by Nicholas Collins on 08/07/2011.
//  Copyright 2011 Nick Collins. All rights reserved.
//
//this code is under GNU GPL3 license, see file COPYING enclosed

#include "PeriodicityAnalysis.h"

#define UKSIZE 256


PeriodicityAnalysis::PeriodicityAnalysis(int maxvoices) {
    
    Ur_ = new float[UKSIZE]; 
    Ud_ = new float[UKSIZE]; 

    maxvoices_ = maxvoices;
    
    voicesaliences_ = new float[maxvoices_]; 
    voiceperiods_ = new float[maxvoices_];
    
    maxq_ = 50; //20 iterations at most; actually, in tests, often needed to go above 20, up to about 40
    
    smax_ = new float[maxq_];
    torlow_= new float[maxq_];
    torup_ = new float[maxq_];
    torprec_ = 0.0000001f; //0.00005f; //start here and increase later
    
    //65 Hz to 2.1KHz allowed search periodicities
    tormin_ = 1.0/2100.0; 
    tormax_ = 1.0/65.0;
    
    K_ = 4096.0/g_samplingrate; 
    
    numvoicesdetected_=0; 
    
    cancellationweight_ = 1.0; 
    polyphonyestimategamma_ = 0.66;
    
//    for (int i=0; i<maxvoices_; ++i) {
//        voiceperiods_[i] = 0.f;  
//        voicesaliences_[i]= 0.f; 
//    }
    
    
}

PeriodicityAnalysis::~PeriodicityAnalysis() {
    
    delete [] Ur_; 
    delete [] Ud_; 
    
    delete [] voicesaliences_;
    delete [] voiceperiods_;
    
    delete [] smax_;
    delete [] torlow_; 
    delete [] torup_;
    
}


////index 15 halfway centres on k in question
//float g_hammingwindownorm[31] = { 0.0056585177222615, 7.0026952423199e-05, 0.0062672829200031, 9.6070554225238e-05, 0.0069081852534025, 0.00014016301897436, 0.007357881966193, 0.00022434355190693, 0.0065604196306213, 0.0004209306170505, 0.0019838514699829, 0.0011244659258033, 0.11559343551383, 0.42817348241183, 0.81822361914331, 1.0, 0.81822361914331, 0.42817348241183, 0.11559343551383, 0.0011244659258033, 0.0019838514699829, 0.0004209306170505, 0.0065604196306213, 0.00022434355190693, 0.007357881966193, 0.00014016301897436, 0.0069081852534025, 9.6070554225238e-05, 0.0062672829200031, 7.0026952423199e-05, 0.0056585177222615 };

//index 4 centre
float g_hammingwindownorm[9] = { 0.0011244659258033, 0.11559343551383, 0.42817348241183, 0.81822361914331, 1.0, 0.81822361914331, 0.42817348241183, 0.11559343551383, 0.0011244659258033 };


//for one frame, iterative reduction 
void PeriodicityAnalysis::compute(float * Uk) {
    
    int i; 
    
    for(i=0; i<UKSIZE; ++i) 
        Ur_[i] = Uk[i]; 
    
    for(i=0; i<UKSIZE; ++i) 
        Ud_[i] = 0.0f; 
    
    //int iteration = 0; 
    
    numvoicesdetected_ = 0; 
    
    //reset everything 
    for (i=0; i<maxvoices_; ++i) {
        voiceperiods_[i] = 0.f;  
        voicesaliences_[i]= 0.f; 
    }
    
    
    prevmixturescore_ = 0.0f; 
    mixturescore_ = 0.0f; 
    
    bool keepgoing = true; 
    
    while (keepgoing) {
        
        //fast  search for min of s
        minSearch(Ur_); 
        
        voicesaliences_[numvoicesdetected_] = bestsalience_;
        voiceperiods_[numvoicesdetected_] = winningtor_;
        
        ++numvoicesdetected_; 
        
        mixturescore_ += bestsalience_; 
        
        float testquantity = mixturescore_/(pow(numvoicesdetected_,polyphonyestimategamma_)); 
        
        if((numvoicesdetected_>=maxvoices_) || (testquantity<=prevmixturescore_)) 
            keepgoing = false; 
        else {
            
            //update Ur, Ud for next round
            
            prevmixturescore_ = testquantity; 
            
            float tor = winningtor_; 
            
            int topm = tor*2745.4833984375; //2745.4833984375/(1.0/tor)
            
            //if(topm>20) topm = 20; 
            
            float srovertor = g_samplingrate/tor;
            //no term for m, since m=1 here for torhat
            
            //numerator free of m precalculated
            float weight = (srovertor + 5.f); // /(srovertor + 320.0f);
            
            for (int m=1; m<topm; ++m) {
     
                int partialk = m*K_/tor +0.5f; //centre for Hamming window reduction
                
                //safety check
                if(partialk<=255) {
                
                    float Urweight = Ur_[partialk]; // * weight; //may be better without this weight multiplier! 
                
                    Urweight *= weight/(m*srovertor + 320.0f);
                    
                //Hamming window added to Ud_
                    
                    int lowk = partialk-4;
                    if(lowk<0) lowk=0; 
                    
                    int highk = partialk+4;
                    if(highk>255) highk=255; 
                    
                    for(int j=lowk; j<=highk; ++j) {
                        
                        int hammingindexnow = j - partialk + 4;
                        
                        float val = g_hammingwindownorm[hammingindexnow] * Urweight; 
      
                        Ud_[j] += val; 
                        
//                        if(m==1) {
//                            
//                            std::cout << "check Ud " << Ud_[j] << " " << Uk[j] << std::endl; 
//                            
//                        }
                    }
                    
                }
                
                    
            }
            
            for(i=0; i<UKSIZE; ++i) {
                
                float diff = Uk[i] - (cancellationweight_*Ud_[i]); 
                
                Ur_[i] = diff>0.f? diff: 0.f; 
                
            }
            
            
        }
    }
      
    if(numvoicesdetected_>0)
        --numvoicesdetected_; //since consistently one too many in tests?
    
}


float PeriodicityAnalysis::smax(int q, float * Ur) {
    
    float tor = 0.5f*(torlow_[q] + torup_[q]);
    float deltator = torup_[q] - torlow_[q];
    
    //2745.4833984375 = 10.7666015625*255 = (44100.0/4096) * 255
    
    
    int topm = tor*2745.4833984375; //2745.4833984375/(1.0/tor)
    
    if(topm>20) topm = 20; 
    
    //break early if go outside first 256 bins (later may allow use of wider band mainFFT for higher k ignoring band wise channel data? 
    
    float salience = 0.0f; 
    float srovertor = g_samplingrate/torup_[q];  //g_samplingrate/tor; 
    
    //float K = 4096.0/g_samplingrate; 
    
    for (int m=1; m<topm; ++m) {
        
        int lowk = m*K_/(tor+0.5*deltator) +0.5f; //dd 0.5 since rounding to nearest integer
        int highk = m*K_/(tor-0.5*deltator) + 0.5f;
        
        //indexing safety check
        if((lowk < UKSIZE) && (highk < UKSIZE)) { 
            
            float maxu = Ur[lowk];
            
            for(int i=lowk+1; i<highk; ++i) {
                
                float nowu = Ur[i];
                
                if(nowu>maxu)
                    maxu = nowu; 
                
            }
            
            float w = 1.0/(m*srovertor + 320.0f);
            
            salience += w* maxu; 
            
        }
        
    }
    
    //was srovertor
    salience *= g_samplingrate/torlow_[q] + 5.f; //numerator multiplier from w doesn't depend on m and taken outside 
    
    //smax_[q] = salience; 
    
    return salience; 
    
}


void PeriodicityAnalysis::minSearch(float * Ur) {
    
    q_ = 0; //Klapuri starts at 1, but we're in C here, not MATLAB
    
    torlow_[0] = tormin_; 
    torup_[0] = tormax_; 
    
    qbest_ = 0; 
    
    while ( ( (torup_[qbest_] - torlow_[qbest_]) > torprec_) && (q_<(maxq_-1))) {
        
        ++q_; 
        
        torlow_[q_] = (torlow_[qbest_] + torup_[qbest_])*0.5f; 
        
        torup_[q_] = torup_[qbest_];
        
        torup_[qbest_] = torlow_[q_]; 
        
        //float s1 = 
        smax_[q_] = smax(q_, Ur);  
        smax_[qbest_] = smax(qbest_, Ur);  
        
        int whichq = 0; 
        float maxval = smax_[0]; 
        
        for (int j=1; j<=q_; ++j) {
            
            float valnow = smax_[j]; 
            
            if(valnow>maxval) {
                
                maxval = valnow; 
                whichq= j; 
                
            }
            
        }
        
        qbest_ = whichq; 
        
    }
    
    //std::cout << "q test" << q_<<std::endl; 
    
    winningtor_ =  (torlow_[qbest_] + torup_[qbest_])*0.5f; 
    bestsalience_ = smax_[qbest_]; 
    
}


