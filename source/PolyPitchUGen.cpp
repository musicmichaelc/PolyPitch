//
//  PolyPitchUGen.cpp
//  Multipitch
//
//  Created by Nicholas Collins on 19/07/2011.
//  Copyright 2011 Nick Collins. All rights reserved.
//
//this code is under GNU GPL3 license, see file COPYING enclosed

#include "PolyPitchUGen.h"

//notes transferred from initial C++ command line project 

//amplitude thresholding, no estimations valid if amp of input too low? 

//compare with full range FFTs on rectified signals? 

//profile - biggest time loss seems to be frontend, particularly sample rate conversion, also filtering in auditory channel and filtering in auditory model bands

//amortise for plug in

//auditory weight front end, so de-emphasise bass? 



PolyPitchUGen::PolyPitchUGen(float samplingRate,int maxvoices, World * world): gamma_(70, 65.f, 5200.f, world), polypitch_(maxvoices) {
    
    pos_ = 0; 
    
    amortisationstage_ = 0; 
    
}

PolyPitchUGen::~PolyPitchUGen() {
    
    
    
}



//AMORTIZED CALCULATIONS 
void PolyPitchUGen::compute(float * input, int numSamples) {
    
    int i,j; 
    
    //update pos and input buffer
    
    float * pointer = inputbuffer_ + pos_; 
    
    //std::cout << amortisationstage_<< " " << numSamples << std::endl;
    
//    for (i=0; i<numSamples; ++i) 
//        printf("i %d  val %f ",i, input[i]); 
//    
//    printf("\n"); 
    
    
    for (i=0; i<numSamples; ++i) 
        pointer[i] = input[i]; 
    
    pos_ +=numSamples;              
    
    if(pos_ >= g_ppwindowsize) {
        
        pos_ = 0; 
        
        //needs to be static during amortisation
        for (i=0; i<g_ppwindowsize; ++i) 
            holdbuffer_[i] = inputbuffer_[i]; 
   
//    for (i=0; i<numSamples; ++i) 
//        printf("i %d  val %f ",i, holdbuffer_[i]); 
//    
//    printf("\n"); 
        
        
        amortisationstage_ = 1;
        
        
        //was <=16, ==17
    } else if ((amortisationstage_>=1) && (amortisationstage_<=20)) {
            
        gamma_.compute(holdbuffer_, g_ppwindowsize,amortisationstage_);
            
        ++amortisationstage_;
        
//    } else if(amortisationstage_==16) {
//   
//        gamma_.compute(holdbuffer_, g_ppwindowsize,amortisationstage_);
//        
    //
//        ++amortisationstage_;
//        
    } else if (amortisationstage_==21) {
        
//        for (j = 0; j<256; ++j) {
//            
//            std::cout << gamma_.Uk_[j];
//            
//            if(j<255) std::cout << " ";
//            
//        }
//        std::cout << std::endl; 

        //gamma_.compute(inputbuffer_, g_ppwindowsize);
        
        
        
        
        
        polypitch_.compute(gamma_.Uk_);
        
        
        
        
        
        //++amortisationstage_;
        amortisationstage_ = 0; 
        
        
//        int voices = polypitch_.numvoicesdetected_; 
//        
//        for  (j = 0; j<voices; ++j) {
//            std::cout << "voice check" << j << " " << (1.0/polypitch_.voiceperiods_[j]) << " " << polypitch_.voicesaliences_[j] << std::endl; 
//        }
               
    } 
    
//    else if(amortisationstage_==17) {
//        
//        int voices = polypitch_.numvoicesdetected_; 
//        
////        cout << i << endl;
//        for  (j = 0; j<voices; ++j) {
////            
////            //cout << i << " " << polypitch.winningtor_ << " " << polypitch.bestsalience_ << endl; 
//            std::cout << j << " " << (1.0/polypitch_.voiceperiods_[j]) << " " << polypitch_.voicesaliences_[j] << std::endl; 
////            
//        }
////        
////        //output Uk values obtained for this frame (256 values)
////        
////        out << i << endl; 
////        
////        for (j = 0; j<256; ++j) {
////            
////            out << gamma.Uk_[j];
////            
////            if(j<255) out << " "; 
////        }
////        
////        out << endl; 
////        
////        out << voices << " "; 
////        
////        for  (int j = 0; j<voices; ++j) {
////            
////            //cout << i << " " << polypitch.winningtor_ << " " << polypitch.bestsalience_ << endl; 
////            out  << (1.0/polypitch.voiceperiods_[j]) << " " << polypitch.voicesaliences_[j]; 
////            
////            if(j<(voices-1)) out << " ";
////            
////        }
////        
////        out << endl;
//        
//        //UPDATE OUTPUT DATA
//        
//        amortisationstage_ = 0;   
//    }
    
}


