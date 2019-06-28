/*
 *  GammatoneFilterBank.h
 *
 *  Created by Nick Collins on 05/06/2009.
 *  Copyright 2009 Nick Collins. All rights reserved.
 *
 */
//this code is under GNU GPL3 license, see file COPYING enclosed

#pragma once 

#include "MultipitchHeader.h"
#include "CalculateFFT.h"

//based on V Hohmann Frequency analysis and synthesis using a Gammatone filterbank Acta Acustica vol 88 (2002): 433--442

//4 first order complex bandpass filters in series

//could build rectification, lowpass, downsampling in as well (use http://www.mega-nerd.com/SRC)? Or just set up in separate wrapper?
//int delay; //may later want delay to keep bands in step
class GammatoneComplexBandpass {
	
public:
	double centrefrequency; 
	double bandwidth; 
	double normalisation; 
	double reala, imaga; 
	//double oldreal[4];
	//double oldimag[4]; 
	double oldreal1,oldreal2,oldreal3,oldreal4;
    double oldimag1,oldimag2,oldimag3,oldimag4;
    
	GammatoneComplexBandpass() {
		
        oldreal1 = 0.0; 
        oldreal2 = 0.0; 
        oldreal3 = 0.0; 
        oldreal4 = 0.0; 
        
        oldimag1=  0.0;
        oldimag2 = 0.0; 
        oldimag3 = 0.0;
        oldimag4 = 0.0;
        
//		for(int i=0; i<4; ++i) {
//			oldreal[i]=0.0;  
//			oldimag[i]=0.0;  
//		}
			
		
	}
	
	~GammatoneComplexBandpass() {
		
	}
	
	void calculate(float * input, float * output, int numSamples); 
	
	
};



//API
//SRC_STATE* src_new (int converter_type, int channels, int *error) ;
//SRC_STATE* src_delete (SRC_STATE *state) ;
//
//int src_process (SRC_STATE *state, SRC_DATA *data) ;
//int src_reset (SRC_STATE *state) ;
//int src_set_ratio (SRC_STATE *state, double new_ratio) ;


extern float g_convertoroutput[g_ppdecimatedwindowsize]; 
extern float g_tempbuffer[g_ppwindowsize]; 


//Klapuri transcription state for channel c at time t
struct AuditoryChannel {
  
    float freq_; 
    float compressionfactor_; 
    double a0_, b1_, a1_; //low pass filter parameters 
    double prev_, previn_;
    
    SRC_STATE* convertorstate_; 
    int convertorerror_; 
    SRC_DATA convertordata_; 
    
    float bandfft_[g_ppdecimatedwindowsize]; //stored output power data for local channel decimated FFT
    //float bandfft_[2048]; //stored output power data for local channel decimated FFT
    
    
    /* Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
     Command line: /www/usr/fisher/helpers/mkfilter -Bu -Lp -o 8 -a 1.1791383220e-01 0.0000000000e+00 -l */
    
//#define NZEROS 8
//#define NPOLES 8
//#define GAIN   1.368365300e+04
//    
//float xv[9], yv[9];

    
    
    
    void setup(float freq) {
        
        freq_ = freq; 
        
        // recursion: tmp = (1-p)*in + p*tmp with output = tmp
        //coefficient: p = (2-cos(x)) - sqrt((2-cos(x))^2 - 1) with x = 2*pi*cutoff/samplerate
        //double propfreq = 2.0*M_PI*freq/g_samplingrate;
//        
//        propfreq = 2.f - cos(propfreq); 
//        double x = propfreq - sqrt(propfreq*propfreq-1); 
//        
        //6dB per oct low pass filter coefficient calculation:
//        double x = exp(-propfreq);
//        a0_ = 1.0-x;
//        b1_ = x;
        
        //void SetLPF(float fCut, float fSampling)
       
            double w = 2.0 * g_samplingrate;
            double norm;
            
            double fCut = freq; 
            fCut *= 2.0 * M_PI;
            norm = 1.0 / (fCut + w);
            b1_ = (w - fCut) * norm;
            a0_ = a1_ = fCut * norm;

        previn_ = 0.0; 
        prev_ = 0.0; 
        
        //sample rate conversion using libsamplerate
        
        //SRC_SINC_MEDIUM_QUALITY
        //SRC_SINC_FASTEST
        //SRC_LINEAR    //leads to double image, aliasing
        //SRC_SINC_FASTEST
        convertorstate_ = src_new (SRC_LINEAR,1,&convertorerror_);
        
        //src_set_ratio(state,0.125);
        
        convertordata_.end_of_input = 0; //more data always available 
        convertordata_.input_frames = g_ppwindowsize;
        convertordata_.output_frames = g_ppdecimatedwindowsize;
        convertordata_.src_ratio = 0.125; 
        convertordata_.data_out = g_convertoroutput; 
        convertordata_.data_in = g_tempbuffer;   

//        for (int i=0; i<9; ++i) {
//            
//            xv[i] = 0.0f;
//            yv[i] = 0.0f;
//        }
            
        
        
    }
    
    void compute (float * target, int numsamples) {
        
        
        //http://www.musicdsp.org/archive.php?classid=3#237

        for (int i=0; i<numsamples; ++i) {
            
            
            //full wave rectify , but leave compression till later
            //double input = target[i]; //
            double input = fabs(target[i]);  
            //HWR
            //double input = target[i];  
            //input = input>0.f?input:0.f; 
            
            
            //double val = a0_*input + b1_*prev_;
            
            double val = input*a0_ + previn_*a1_ + prev_*b1_;
            
            previn_ = input; 
            
            target[i] = val;
            
            //original HWR effect pre downsampling and FFT
            //target[i] = 0.5*(val+target[i]);
            
            prev_ = val; 
            
        }

        src_process(convertorstate_,&convertordata_); 
        
        
        //downsampling using low pass filter and drop 7 out of 8
        //Butterworth low pass filters designed via http://www-users.cs.york.ac.uk/~fisher/mkfilter/trad.html
        //CPU cost more expensive than libsamplerate use
        
//        for (int i=0; i<numsamples; ++i) {
//            
////        xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4]; xv[4] = xv[5]; xv[5] = xv[6]; xv[6] = xv[7]; xv[7] = xv[8]; 
////            xv[8] = target[i] *7.3079900520716e-05f; // / 1.368365300e+04;
////        yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; yv[4] = yv[5]; yv[5] = yv[6]; yv[6] = yv[7]; yv[7] = yv[8]; 
////        yv[8] =   (xv[0] + xv[8]) + 8 * (xv[1] + xv[7]) + 28 * (xv[2] + xv[6])
////        + 56 * (xv[3] + xv[5]) + 70 * xv[4]
////        + ( -0.0199403269 * yv[0]) + (  0.2350845773 * yv[1])
////        + ( -1.2374188036 * yv[2]) + (  3.8086314687 * yv[3])
////        + ( -7.5233087140 * yv[4]) + (  9.8112487180 * yv[5])
////        + ( -8.3035924287 * yv[6]) + (  4.2105870548 * yv[7]);
////        
//            
//            //4th order
//            xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4]; 
//            xv[4] = target[i] * 0.0084081014081739f;
//            yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; 
//            yv[4] =   (xv[0] + xv[4]) + 4 * (xv[1] + xv[3]) + 6 * xv[2]
//            + ( -0.1366748840 * yv[0]) + (  0.8077385823 * yv[1])
//            + ( -1.8873902282 * yv[2]) + (  2.0817969073 * yv[3]);
//            
//            //assumes numsamples is divisible by 8
//           if(i%8==0)  
//           g_convertoroutput[i/8] = yv[4]; //yv[8]
//            
//        }
        
        //results ready for FFT in g_convertoroutput
        
    }
    
    
};


class GammatoneFilterBank {
	
public:
	int numbands; 
	GammatoneComplexBandpass * bands; 
    
	//could be global of whole system
	double samplingperiod; //T 
	
    AuditoryChannel * channels_;
    //float * mainfft_; 
    //float * decimatedfft_; 
    //float overallscalefactor; 
    
    float * Uk_; //Klapuri notation  
    int * frprevchannel_; 
    int * frchannel_; 
    float * frinterp_; 
    
    float levelcompressionfactor_; 
    float mixleftterm_;
    
    CalculateFFT * mainfft_; 
    CalculateFFT * decimatedfft_; 
    
    float mainfftoutput_[g_ppwindowsize];
    
	GammatoneFilterBank(int numbands=88, float lowfreq=100.0, float highfreq=18000.0, World * world=0): numbands(numbands) {
		
        levelcompressionfactor_ = -0.1; 
        mixleftterm_ = 4.0; 
        
		bands= new GammatoneComplexBandpass[numbands];
		
		//(samplingrate)
			
		samplingperiod= g_sampleperiod; //1.0/samplingrate;
		
		designFilterBank(lowfreq, highfreq); 
		
        channels_ = new AuditoryChannel[numbands]; 
        
        for (int i=0; i<numbands; ++i) 
            channels_[i].setup(bands[i].centrefrequency); 
		  
        //zero padding, Hann window
        mainfft_ = new CalculateFFT(g_ppwindowsize,2*g_ppwindowsize,0,world);
        
        //decimatedfft_ = new CalculateFFT(2048,4096); //old supposedly "less" efficient version, full FFT on each band
        decimatedfft_ = new CalculateFFT(g_ppdecimatedwindowsize,2*g_ppdecimatedwindowsize,0,world);
        
        //final periodicity ratings from 65 to 2.1 kHz via 2048 FFT at R=44100Hz
        Uk_ = new float[g_Uksize]; //same size as zero padded FFT
        
        //frequencyresponse_ = new float[256]; 
        
        frprevchannel_ = new int[g_Uksize];
        frchannel_ = new int[g_Uksize];
        frinterp_ = new float[g_Uksize];
        initFrequencyResponse(); 
        
	}
	
	~GammatoneFilterBank() {
		
		delete [] bands; 
		
        delete [] channels_; 
        
        delete mainfft_; 
        
        delete decimatedfft_; 
        
        delete [] Uk_; 
        delete [] frprevchannel_; 
        delete [] frchannel_; 
        delete [] frinterp_; 
	}
	
    void initFrequencyResponse(); 
    
	void designFilterBank(float lowfreq, float highfreq); 
	
    void compute(float * input, int numsamples, int amortisationstage);
    
};


