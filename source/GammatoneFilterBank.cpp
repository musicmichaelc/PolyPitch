/*
 *  GammatoneFilterBank.cpp
 *  Created by Nick Collins on 05/06/2009.
 *  Copyright 2009 Nick Collins. All rights reserved.
 *
 */

//this code is under GNU GPL3 license, see file COPYING enclosed


//IEM outer and middle ear filters
//http://iem.at/projekte/publications/iem_report/report02_98/
//some generally good stuff there! 
//also see PH96 

#include "GammatoneFilterBank.h"

#include <math.h> 
#include <stdio.h> 

float g_convertoroutput[g_ppdecimatedwindowsize]; 
float g_tempbuffer[g_ppwindowsize]; 

//float g_mainfftoutput[2048]; //half fft size


//float g_tempbuffer2[2048]; 

//bandwise updates, in window of numsamples 

//will probably want a way to amortise over bands within a frame: 64 lots of 64 within 4096, but could do in batches
//of X bands at a time, calculated from knowledge of working blocksize
void GammatoneFilterBank::compute(float * input, int numsamples, int amortisationstage) {
    
    int i,j; 
	
    float * output; 
    
  
    if (amortisationstage==1) {
   
    mainfft_->calculateFrame(input,mainfftoutput_); 
    
        //was 16
    } else if (amortisationstage<20){
 
        //70 bands over 14 amortisationstages, 5 per stage
        //int startj = (amortisationstage-2)*5;
        //int endj = startj+5;
        
        //still too heavy peak load, so 
        //70 bands over 18 amortisationstages, 4 per stage except 18th just 2
        
        int startj = (amortisationstage-2)*4;
        int endj = startj+4;
        
        if(amortisationstage==19)
           endj = startj+2; 
        
        
        //for (j=0; j<numbands; ++j) 
        for (j=startj; j<endj; ++j) {
            
            bands[j].calculate(input, g_tempbuffer, numsamples); 
            
            //do something further with each band
            output = g_tempbuffer; 
            
            //compression: find standard deviation in this frame
            
            float meansum = 0.0f; 
            float squaresum = 0.0f; 
            //find mean
            
            //and run cheaper standard dev equation E(X^2)- mean^2 simultaneously?
            
            float valuenow;
            
            for (i=0; i<numsamples; ++i) {
                
                valuenow = output[i]; 
                
                meansum += valuenow; 
                squaresum += valuenow*valuenow; 
                
            }
            
            meansum /= numsamples; 
            
            float variance = squaresum/numsamples - meansum*meansum; 
            
            //combine stddev = sqrt(variance) with Klapuri equation (5)
            //smaller factors actually work better, esp. at preserving differences of amplitude in original rather than over compressing
            
            float compression = powf(variance,levelcompressionfactor_); //powf(variance,-0.1); //1.0f; //powf(variance,-0.33); //-0.33 = 0.5*(-0.66) = 0.5* (0.33-1)
            
            //apply compression and half wave rectify simultaneously to avoid two loops through
            //        for (i=0; i<numsamples; ++i) {
            //            
            //            valuenow = output[i]; 
            //            
            //            //HWR
            //            //output[i] = valuenow>0.f?(valuenow*compression):0.f; 
            //            //full wave rectification
            //            
            //            valuenow *=compression; 
            //            
            //            output[i] = valuenow;
            //            
            //            g_tempbuffer2[i] = fabs(valuenow); 
            //             
            //        }
            
            //should low pass filter below 2cf
            
            //Klapuri uses a different sequence of steps, p258 end of section B
            
            //isn't low pass filter equivalent to taking FFT and discounting certain frequency components before 
            //doing power spectrum, IFFT? Perhaps not because of zero padding? Hmmm
            
            //could do brute force comb/time domain auto corr only at piano note locations/lags? 
            
            
            channels_[j].compressionfactor_ = compression; 
            
            
            //printf("compression %d %f \n", j, channels_[j].compressionfactor_); 
            
            //update lpf and do sample rate conversion
            //low pass filter and decimate by factor of 8, FFT every 256 samples (main FFT every 2048)
            
            channels_[j].compute(output, numsamples); 
            //channels_[j].compute(g_tempbuffer2, numsamples);
            
            decimatedfft_->calculateFrame(g_convertoroutput,channels_[j].bandfft_); 
            
            
            //less efficient way (apparently...)
            //        for (i=0; i<numsamples; ++i) 
            //            output[i] = 0.5*(output[i] + g_tempbuffer2[i]); 
            //        
            //        decimatedfft_->calculateFrame(output,channels_[j].bandfft_); 
            
        }
        

        
        
    } else {
        
        //now calculate Uk over 256 bands
        //factor of 2 dropped, p=1, so just fabs; are using powers for x, though... so may need to restore sqrt (in CalculateFFT probably) 
        
        //lefttermmult = 1.0f to start  with? 
        //float lefttermmult = 4.0f; //4.0f
        
        //    int channelnow = 0; 
        //    float channelfreq; 
        //    int topchannel = numbands-1; 
        //    int channelprev; 
        
        //float hzperbin = g_sampleperiod/4096.f; 
        
        //should precalculate frequency response rather than each time? 
        
        //initFrequencyResponse(); //actually needs recalculating each time, bollocks! 
        //could store interpolation positions? 
        
        float prevcomp; 
        float nextcomp; 
        float t, interp; 
        
        for (i=1; i<256; ++i) {
            
            // float freqnow = i*hzperbin;  
            
            //sum up channel FFTs
            
            float sum = 0.0f; 
            
            for (j=0; j<numbands; ++j) {
                
                //sum += fabs(channels_[j].compressionfactor_* channels_[j].bandfft_[i]);
                
                sum += channels_[j].compressionfactor_* channels_[j].bandfft_[i];
                
                //sum += channels_[j].bandfft_[i];
            }
            
            //printf("i %d sum bands %f ",i, sum); 
            
            prevcomp = channels_[frprevchannel_[i]].compressionfactor_;
            nextcomp = channels_[frchannel_[i]].compressionfactor_;
            t  =  frinterp_[i];  
            interp = (1.0-t)*prevcomp + t*nextcomp; 
            
            //lefttermmult
            sum += mixleftterm_* interp * mainfftoutput_[i] ; //fabs(interp * g_mainfftoutput[i]); 
            
            //printf("sum main %f ",sum, mainfftoutput_[i]); 
            
            //if(i<6) sum = 0.0f; //no ability to cope with f0 under 65Hz, around bin 6
            
            Uk_[i] = sum; 
        }
        
        Uk_[0]= 0.0f; 
        
    }

    
    
  
    
}



void GammatoneFilterBank::initFrequencyResponse() {
    
    int i; 
    float hzperbin = g_samplingrate/4096.f; 
    
    int channelnow = 0; 
    float channelfreq; 
    float channelprevfreq; 
    int topchannel = numbands-1; 
    int channelprev; 
    float interp; 
    
    for (i=1; i<256; ++i) {
        
        float freqnow = i*hzperbin;  
        
        //update bandbelow to be just below or equal to freqnow
        
        bool keepgoing = true;
        
        while(keepgoing) {
            
            channelfreq = channels_[channelnow].freq_; 
            
            if((channelfreq>freqnow) || (channelnow==topchannel)) {
                keepgoing = false;
                
            } else
                ++channelnow; 
            
        }
        
        //interpolation of compressionfactor_ function
        
        //float compressioninterp = 1.0f; 
        
        if(channelnow==0) {
            //compressioninterp = channels_[0].compressionfactor_; 
            channelprev=0;
            interp = 0.0; 
        
        }
        else {
            
            if (freqnow>channelfreq) {
                //compressioninterp = channels_[topchannel].compressionfactor_; 
                channelprev = topchannel; 
                interp = 0.0; 
                
            }
                else {
                   
                    //proper inbetween interpolation
                    
                    
                    channelprev = channelnow-1; //valid since channelnow at least 1
                    
                    channelprevfreq = channels_[channelprev].freq_;
                    
                    //linear, oh well, will do
                    interp = (freqnow - channelprevfreq)/(channelfreq - channelprevfreq); 
                    
                }
                
            
        }
        
        frprevchannel_[i] =  channelprev;
        frchannel_[i] =  channelnow;
        frinterp_[i] =  interp;  
        //frequencyresponse_[i] = compressioninterp; 
        
    }
}




void GammatoneFilterBank::designFilterBank(float lowfreq=50.0, float highfreq=20000.0) {
	
	int i; 
	
	//1 erb= 26.014253045474Hz
	
	//from Klapuri 2006
	//f = 229*(10**(0.046728972*e)-1)
	//e= log10((f/229 +1))/0.046728972
	
	double lowerbs = log10((lowfreq*0.00437 +1))*21.4; 
	double higherbs = log10((highfreq*0.00437 +1))*21.4; 
	double prop= 0.0;
	double erbsnow; 
	double cf, b; 
	double beta, phi, p, lambda; 
		
	for(i=0; i<numbands; ++i) {
		
		prop= ((double)i)/(numbands-1);
		
		erbsnow= lowerbs+((prop)*(higherbs-lowerbs)); 
		
		GammatoneComplexBandpass * band= &(bands[i]); 
		
		cf=  228.8455*(pow(10, 0.046728972*erbsnow)-1);
		
		//get cfs from externally calculated 12TET data for 88 piano notes 
		
		//cf= g_pianonotecfs[i]; 
		
		
		//*2pi
		//b= 158.143375952
		//2*pi*1.019*24.7*(4.37*cf/1000+1)
		//b in Hz not angular frequency	for calculations in Hohmann
		b= 24.7*(cf*0.00437 +1);
		
		//reducing bandwidth as test of resonance 
		//b *=0.125; 
		//printf("band num %d erbsnow %f lowerbs %f higherbs %f cf %f b %f\n", i, erbsnow, lowerbs, higherbs, cf, b); 
		
		band->centrefrequency =cf; 
		
		//actually need to convert ERBs to 3dB bandwidth
		b= 0.887*b; //converting to 3dB bandwith in Hz, 	//PH96 pg 3
		
		band->bandwidth= b; 
	
		// filter coefficients to calculate, p.435 hohmann paper
		
		beta= 6.2831853071796*cf*samplingperiod;
		phi= 3.1415926535898*b*samplingperiod; 
		p=  (1.6827902832904*cos(phi) -2)*6.3049771007832;  
		lambda= (p*(-0.5))-(sqrt(p*p*0.25-1.0)); 
		
		band->reala= lambda*cos(beta); 
		band->imaga= lambda*sin(beta);
		
		//avoid b= 0 or Nyquist, otherise must remove factor of 2.0 here
		band->normalisation= 2.0*(pow(1-fabs(lambda),4)); 
		
		//printf("band %d norm %.12f \n", i, band->normalisation); 
		//works perfectly according to equations in paper; so search for first max? 
		//if(i>2 && i<10) {
			
	//		double norm= band->normalisation; 
//			
//			//printf("impulse response %d lambda %f powertest %f norm %.12f \n", i, lambda, pow(lambda,50), norm); 
//			
//			for (j=0; j<1000; ++j)
//			g_debugcheck[j]=0.0; 
//			
//			g_debugcheck[0]=1.0; 
//			
//			band->calculate(g_debugcheck,g_debugcheck,1000); 
//			
//			//manually find envelope maximum
//			float maxima=0.0; 
//			int maxind=0; 
//				
//			for (j=0; j<1000; ++j) {
//			
//				float nextval= fabs(g_debugcheck[j]); 
//				
//				if(nextval>maxima) {maxima= nextval; maxind= j;}
//				
//			}
//			
//			//
////			//printing impulse response
////			for (j=0; j<1000; ++j) {
////				
////				float poly= ((j*j*j)+ (6*j*j) + (11*j) +6); 
////				float real= pow(lambda,j)*cos(beta*j);
////				
////				printf("m: %d theoretical: %f actual: %f \t", j, norm*real*poly/6.0, g_debugcheck[j]); 
////				
////			}
//			
//			printf("max of impulse response for band i %d =  %d \n", i, maxind); 
//			
			
		//}
		
	}
	
}




// this is a function for preventing pathological math operations in ugens.
// can be used at the end of a block to fix any recirculating filter values.
//inline float zapgremlins(float x)
//{
//	float absx = fabs(x);
//	// very small numbers fail the first test, eliminating denormalized numbers
//	//    (zero also fails the first test, but that is OK since it returns zero.)
//	// very large numbers fail the second test, eliminating infinities
//	// Not-a-Numbers fail both tests and are eliminated.
//	return (absx > (float)1e-15 && absx < (float)1e15) ? x : (float)0.;
//}


//can at least save one add in imag part of first filter due to real signal input
//optimisations; just declare all as single variables rather than dereference pointers? 
void GammatoneComplexBandpass::calculate(float * input, float * output, int numSamples) {
	
	int i,j; 
	
	double newreal, newimag; 
	
	
    //loop unwinding, avoid array access
    
    
    for (i=0; i<numSamples; ++i) {
        
		newreal= input[i]; //real input 
		newimag=0.0; 
		
        newreal += (reala*oldreal1)-(imaga*oldimag1);
        newimag += (reala*oldimag1)+(imaga*oldreal1);
        
        oldreal1= newreal; 
        oldimag1= newimag; 
        
        newreal += (reala*oldreal2)-(imaga*oldimag2);
        newimag += (reala*oldimag2)+(imaga*oldreal2);
        
        oldreal2= newreal; 
        oldimag2= newimag; 
        
        newreal += (reala*oldreal3)-(imaga*oldimag3);
        newimag += (reala*oldimag3)+(imaga*oldreal3);
        
        oldreal3= newreal; 
        oldimag3= newimag; 
        
        newreal += (reala*oldreal4)-(imaga*oldimag4);
        newimag += (reala*oldimag4)+(imaga*oldreal4);
        
        oldreal4= newreal; 
        oldimag4= newimag; 
        
        
//		for (j=0; j<4; ++j) {
//			
//			newreal= newreal + (reala*oldreal[j])-(imaga*oldimag[j]);
//			newimag= newimag + (reala*oldimag[j])+(imaga*oldreal[j]);
//			
//			oldreal[j]= newreal; //zapgremlins(newreal); //trying to avoid denormals which mess up processing via underflow
//			oldimag[j]= newimag; //zapgremlins(newimag); 
//		}
		
		output[i]= newreal*normalisation; 
		
		//imaginary output too could be useful
		
	}

    
//old version
    //float newreal, newimag; 
//	for (i=0; i<numSamples; ++i) {
//	
//		newreal= input[i]; //real input 
//		newimag=0.0; 
//		
//		for (j=0; j<4; ++j) {
//			
//			newreal= newreal + (reala*oldreal[j])-(imaga*oldimag[j]);
//			newimag= newimag + (reala*oldimag[j])+(imaga*oldreal[j]);
//			
//			oldreal[j]= newreal; //zapgremlins(newreal); //trying to avoid denormals which mess up processing via underflow
//			oldimag[j]= newimag; //zapgremlins(newimag); 
//		}
//		
//		output[i]= newreal*normalisation; 
//		
//		//imaginary output too could be useful
//		
//	}
	
}













//SC code
//Array.fill(88, {|i| (21+i).midicps}).postcs; 
//
//(43..96).size  //54
//
//Array.fill(54, {|i| (43+i).midicps}).postcs; 
//float g_pianonotecfs[88]= {27.5, 29.135235094881, 30.867706328508, 32.703195662575, 34.647828872109, 36.708095989676, 38.89087296526, 41.203444614109, 43.653528929125, 46.249302838954, 48.999429497719, 51.913087197493, 55.0, 58.270470189761, 61.735412657016, 65.40639132515, 69.295657744218, 73.416191979352, 77.78174593052, 82.406889228217, 87.307057858251, 92.498605677909, 97.998858995437, 103.82617439499, 110.0, 116.54094037952, 123.47082531403, 130.8127826503, 138.59131548844, 146.8323839587, 155.56349186104, 164.81377845643, 174.6141157165, 184.99721135582, 195.99771799087, 207.65234878997, 220.0, 233.08188075904, 246.94165062806, 261.6255653006, 277.18263097687, 293.66476791741, 311.12698372208, 329.62755691287, 349.228231433, 369.99442271163, 391.99543598175, 415.30469757995, 440.0, 466.16376151809, 493.88330125612, 523.2511306012, 554.36526195374, 587.32953583482, 622.25396744416, 659.25511382574, 698.45646286601, 739.98884542327, 783.9908719635, 830.60939515989, 880.0, 932.32752303618, 987.76660251225, 1046.5022612024, 1108.7305239075, 1174.6590716696, 1244.5079348883, 1318.5102276515, 1396.912925732, 1479.9776908465, 1567.981743927, 1661.2187903198, 1760.0, 1864.6550460724, 1975.5332050245, 2093.0045224048, 2217.461047815, 2349.3181433393, 2489.0158697766, 2637.020455303, 2793.825851464, 2959.9553816931, 3135.963487854, 3322.4375806396, 3520.0, 3729.3100921447, 3951.066410049, 4186.0090448096};
//
////MIDI 43 to 96 inclusive
//float g_otherbands[54]= { 97.998858995437, 103.82617439499, 110.0, 116.54094037952, 123.47082531403, 130.8127826503, 138.59131548844, 146.8323839587, 155.56349186104, 164.81377845643, 174.6141157165, 184.99721135582, 195.99771799087, 207.65234878997, 220.0, 233.08188075904, 246.94165062806, 261.6255653006, 277.18263097687, 293.66476791741, 311.12698372208, 329.62755691287, 349.228231433, 369.99442271163, 391.99543598175, 415.30469757995, 440.0, 466.16376151809, 493.88330125612, 523.2511306012, 554.36526195374, 587.32953583482, 622.25396744416, 659.25511382574, 698.45646286601, 739.98884542327, 783.9908719635, 830.60939515989, 880.0, 932.32752303618, 987.76660251225, 1046.5022612024, 1108.7305239075, 1174.6590716696, 1244.5079348883, 1318.5102276515, 1396.912925732, 1479.9776908465, 1567.981743927, 1661.2187903198, 1760.0, 1864.6550460724, 1975.5332050245, 2093.0045224048}; 

//float g_debugcheck[1000]; 




