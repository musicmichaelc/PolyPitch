/*
 *  CalculateFFT.h
 *
 *  Created by Nick Collins on 23/11/2008.
 *  Copyright 2008 Nick Collins. All rights reserved.
 *
 */

//this code is under GNU GPL3 license, see file COPYING enclosed

#pragma once


#include "MultipitchHeader.h"

#ifdef OSXACCELERATENOTSCFFT
#include <Accelerate/Accelerate.h>
#else
#include "SC_fftlib.h"
#endif


//#include "FFT_UGens.h"

//#include <vDSP.h>
//#include "vecLib/vDSP.h"

//does this deal with its own hop and shunt? 



class CalculateFFT {
    
private:
    float * fftbuffer;
	float * zeropad; 	
    //float * prefftbuffer;
    //float * postfftbuffer;
    float * window;
    
    
    
#ifdef OSXACCELERATENOTSCFFT
    //vDSP
    unsigned long vlog2n;
    COMPLEX_SPLIT splitBuf;
    FFTSetup vsetup;
#else    
    scfft* m_scfft;
#endif
    
    void prepareHannWindow();
    void prepareHammingWindow();
    
public:
    //default to hanning window for now
    //CalculateFFT(int fftsize, int hopsize);
    //~CalculateFFT();
    
    int fftsize;
    int inputsize; 
    int padding; 
    int nover2;
    int dowindowing; 
    float reciprocalfftsize; 
    
    World * world_;     
    
    //int hopsize;
    //int halfsize;
    
    //could have a Spectrum object to hold this data, and override for other spectral representations? 
    //returns power spectrum
    void calculateFrame(float * input, float * output);
    
    //CalculateFFT::CalculateFFT(int fftsize, int hopsize) : fftsize(fftsize), hopsize(hopsize) {
    //int windowflag=0, dowindowing(windowflag)
    CalculateFFT(int inputsize=1470, int fftsize=2048, float * externalwindow=0, World * world=0) : inputsize(inputsize), fftsize(fftsize), window(externalwindow) {
        
        
        //printf("CalculateFFT %p \n",world); 
        
        nover2= fftsize >>1;
        
        reciprocalfftsize = 1.0f/((float)fftsize); 
        
        //in place
        //fftbuffer= new float[fftsize];
        
        //prefftbuffer = new float[fftsize];
        //postfftbuffer = new float[fftsize];
        
        //externalwindow
        
        padding = fftsize-inputsize; 
        //	
        //
        //	//not needed
        //	zeropad= new float[fftsize];
        //
        //	for (int i = 0; i < fftsize; ++i) {
        //		
        //		zeropad[i]=0.0; //assuming 2048-1470 end samples of this should stay zero always
        //	
        //	}
        //	
        //	
        if(externalwindow ==0) {
            
            window= new float[inputsize];
            
            prepareHammingWindow(); 
        }
        
        //prepareHanningWindow();
        
        
        fftbuffer = ( float* ) malloc ( fftsize * sizeof ( float ) ); 
        
#ifdef OSXACCELERATENOTSCFFT
        
        ////////vDSP///////////////
        // Allocate memory for the input operands and check its availability, 
        // use the vector version to get 16-byte alignment; required for vectorized code 
        //used to be vec_malloc, deprecated, now malloc
        
        splitBuf.realp = ( float* ) malloc ( nover2 * sizeof ( float ) ); 
        splitBuf.imagp = ( float* ) malloc ( nover2 * sizeof ( float ) ); 
        
        int vlog= 0; 
        
        int currfftsize= fftsize; 
        // printf("currfftsize %d vlog %d \n", currfftsize, vlog);
        while(currfftsize>1) {++vlog; currfftsize>>=1;}
        //printf("log of n =%d is %d \n", fftsize, vlog); 
        
        vlog2n = vlog; //11; //hard code for now log2(fftsize); //10; //N is hard coded as 1024, so 10^2=1024 //log2max(N);
        //vsetup = create_fftsetup(vlog2n, 0);
        
        vsetup= vDSP_create_fftsetup(vlog2n, FFT_RADIX2);
        
        //	planTime2FFT = fftwf_plan_r2r_1d(fftsize, prefftbuffer, postfftbuffer, FFTW_R2HC, FFTW_ESTIMATE);	
        
        //halfsize= fftsize/2;
        //obtainedReal = ( float* ) malloc ( fftsize * sizeof ( float ) ); 
        
        //splitBuf.realp = new float [nover2]; //  * sizeof(float)); 
        //splitBuf.imagp = new float [nover2]; //(float*)RTAlloc(unit->mWorld, nover2 * sizeof(float));
        
        //get LOG2CEIL as per sc_fft
        
        
        
        
#else
        
        //printf("CalculateFFT 2 %p \n",world); 
        
        SCWorld_Allocator alloc(ft, world);
        
        //printf("CalculateFFT 3a %p %p %p %d %d %p \n",ft, world, fftbuffer, fftsize, inputsize, &alloc); 
        
        //rectangular, will pre-multiply with hamming ourselves and do zero padding
        m_scfft = scfft_create(fftsize, fftsize, kRectWindow,fftbuffer,fftbuffer, kForward, alloc);
        
        //printf("CalculateFFT 3b %p \n",world); 
        
#endif	
        
        world_ = world; 
        
        
    }
    
    
    ~CalculateFFT() {
        
#ifdef OSXACCELERATENOTSCFFT
        
        if (vsetup) vDSP_destroy_fftsetup(vsetup);
        if (splitBuf.realp) free(splitBuf.realp); //delete [] splitBuf.realp;
        if (splitBuf.imagp) free(splitBuf.imagp); //delete [] splitBuf.imagp;
        
#else
        SCWorld_Allocator alloc(ft, world_);
        
        if(m_scfft) 
            scfft_destroy(m_scfft, alloc);
#endif	   
        
        free(fftbuffer);
        
        //delete [] fftbuffer;
        //delete [] prefftbuffer;
        ///delete [] postfftbuffer;
        delete [] window;
        //delete [] zeropad; 
        
    }
    
    
    
    
};