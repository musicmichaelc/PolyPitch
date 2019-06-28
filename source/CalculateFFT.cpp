/*
 *  CalculateFFT.cpp
 *  openFrameworks
 *
 *  Created by Nick Collins on 23/11/2008.
 *  Copyright 2008 Nick Collins. All rights reserved.
 *
 */

//this code is under GNU GPL3 license, see file COPYING enclosed

#include "CalculateFFT.h"
#include <math.h>
#include <iostream>

#define TWOPI 6.28318530717952646f 



//probably have to put this in load() again since 256 or more cos calls
void CalculateFFT::prepareHannWindow() {
	float ang;
	
	ang=(1.0/inputsize)*TWOPI; //2.0*M_PI;
	
	for(int i=0;i<inputsize;++i)
		window[i]=0.5 - 0.5*cos(ang*i);
	
}

void CalculateFFT::prepareHammingWindow() {
	float ang;
	
	ang=(1.0/(inputsize-1))*TWOPI; //2.0*M_PI;
	
	for(int i=0;i<inputsize;++i)
		window[i]=0.54 - 0.46*cos(ang*i);
	
}

//hopefully no problems with max size of int 
//INT_MAX in limits.h 2147483648/44100/60 = 811.59623847317 minutes of audio in a file

//does all at once, could convert to just one frame at a time? 
void CalculateFFT::calculateFrame(float * input, float * output) {
	
	int j; //i, pos;
	
	//copy first	
	if(window!=0) {
        for (j=0; j<inputsize; ++j)
			fftbuffer[j] = input[j]* window[j];
	} else {
        for (j=0; j<inputsize; ++j)
            fftbuffer[j] = input[j];  
	}
    
	for (j=0; j<padding; ++j)
		fftbuffer[inputsize+j] = 0.0f;  
	
    
#ifdef OSXACCELERATENOTSCFFT    
    
	//go back to vDSP instructions? 
	//look for most straight forward real FFT using vDSP
	
	
	//stages in vDSP for real vector: 
	//1) transform real input to odd then even vector (ctoz)
	//2) do complex fft
	//3) return to real format (ztoc)
	
    //	// Look at the real signal as an interleaved complex vector by casting it.
    //    // Then call the transformation function ctoz to get a split complex vector,
    //    // which for a real signal, divides into an even-odd configuration.
    //    ctoz ((COMPLEX *) fftbuf, 2, &unit->m_vA, 1, NOVER2);
    //	
    //    // Carry out a Forward FFT transform
    //    fft_zrip(unit->m_vsetup, &unit->m_vA, 1, unit->m_vlog2n, FFT_FORWARD);
    //	
    //    // The output signal is now in a split real form.  Use the function
    //    // ztoc to get a split real vector.
    //    ztoc ( &unit->m_vA, 1, (COMPLEX *) fftbuf, 2, NOVER2);
    //		
    
	// Perform even-odd split
	vDSP_ctoz((COMPLEX*) fftbuffer, 2, &splitBuf, 1, nover2);
	// Now the actual complex FFT; origin of extra scale factor of 2 is this! 
	vDSP_fft_zrip(vsetup, &splitBuf, 1, vlog2n, FFT_FORWARD);
	// Copy the data to the public output buf, transforming it back out of "split" representation
	vDSP_ztoc(&splitBuf, 1, (COMPLEX*) fftbuffer, 2, nover2);
	//was in place
    
	//undo scale factor of 2 caused by complex fft on 2n 
	//can copy to output buffer at the same time? 
	
	float scale= 0.5f/fftsize;  //compensate for FFT size in scale factor
    
	//vDSP_vsmul( fftbuffer, 1, &scale, output, 1, fftsize); 
	vDSP_vsmul( fftbuffer, 1, &scale, fftbuffer, 1, fftsize); 
    
#else
    
    scfft_dofft(m_scfft);
    
    //scaling
    for(j=0; j<fftsize; ++j) 
        fftbuffer[j] *= reciprocalfftsize; 
    
#endif
	
	output[0]= sqrt(fftbuffer[0]* fftbuffer[0]); //get power
    //	output[0]= 10*log10((953674*output[0])+1);
    //	
    
    int halfsize = fftsize/2; 
    
    
	// Squared Absolute so get power
	//do half of the calculations
	for(j=1; j<halfsize; j++) {
        int index = 2*j; 
		float val1= fftbuffer[index];
		float val2= fftbuffer[index+1];
		
		output[j] = sqrt((val1*val1) + (val2*val2)); //magnitudes, not powers 
		
		//output[j]= 10*log10((953674*output[j])+1);
		
		//should do log conversions as 0 to 120 dB (16 bit range- but some samples are 24 bit!) where x is a power after fft, 
		//multiplier for 96dB range = 3797 
		//120 dB range = 953674
		//logpower= 10*log10(3797*x+1)
		
		//log conversions up to 105 dB max  (minimum is -30 however!)
		//float db= 10*log10((bsum*32382)+0.001); 
		
		//if(j<5) printf("j %d val %f \t",j, output[j]);
		//fftbuf[256-j] = 0.0f;
	}
	
	
	
}
//elements 0 to 512 of array2 first 513 bins/spectral buckets, at 0, 43, 86, ..., 22050 Hz












