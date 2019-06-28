//SuperCollider is under GNU GPL version 3, http://supercollider.sourceforge.net/
//these extensions released under the same license

/*
 *  PolyPitch.cpp
 *
 *  Created by Nicholas Collins on 19/07/2011
 *  Copyright 2011 Nicholas M Collins. All rights reserved.
 *
 */
//this code is under GNU GPL3 license, see file COPYING enclosed


//building with CMake
//need additional libs in both architectures if linking in 
//cmake -DSC_PATH=/data/gitprojects/supercollider -DCMAKE_OSX_ARCHITECTURES='x86_64' ..

//todo:
//test SinOsc inputs picked up on...
//back to C code with some clearly periodic test signals



//#include "MultipitchHeader.h"

#include "PolyPitchUGen.h"

InterfaceTable *ft; 

float g_samplingrate = 44100.f; 
float g_sampleperiod = 2.2675736961451e-05f; 


//parameters: 
//maxq_, torprec_, lefttermmult, -0.1 compression constant

struct PolyPitch : public Unit  
{
    PolyPitchUGen * ugen_;
    int m_blocksize; 
    int m_sr; 
    
    int maxvoices_; 
    float levelcompressionfactor_; 
    float mixleftterm_; 
    float torprec_; 
    float cancellationweight_; 
    float polyphonyestimategamma_;
    
};



//data to be shared between RT and NRT threads 
struct CmdData {
	enum Type
    {
		NRTCtorPP,
		NRTDtorPP
    };
	
	Type type; 
	Unit * unit; 
	void *	nrtallocated; 
	float samplingrate_; 
	
};


extern "C" {  
	
	void PolyPitch_next(PolyPitch *unit, int inNumSamples);
	void PolyPitch_Ctor(PolyPitch* unit);
	void PolyPitch_Dtor(PolyPitch* unit);


}


//since classes are simple, can use basic placement new and explicit destructor call 

//for NRT allocation and deallocation 

bool cmdStage2(World* inWorld, CmdData* cmd) // NRT
{
	//Unit* unit = cmd->unit;
	
	switch (cmd->type) {
		case CmdData::NRTCtorPP: {
            
            //set up algorithm values
            PolyPitch* pUnit = (PolyPitch*)cmd->unit;
   
            PolyPitchUGen * pPoly = new PolyPitchUGen(cmd->samplingrate_,pUnit->maxvoices_,inWorld); 
            
            pPoly->gamma_.levelcompressionfactor_ = pUnit->levelcompressionfactor_;
            pPoly->gamma_.mixleftterm_ = pUnit->mixleftterm_;
            
            pPoly->polypitch_.torprec_ = pUnit->torprec_;
            pPoly->polypitch_.cancellationweight_ = pUnit->cancellationweight_;
            pPoly->polypitch_.polyphonyestimategamma_ = pUnit->polyphonyestimategamma_;
 
			cmd->nrtallocated = (void *)pPoly;
            
        }
			return true;
		case CmdData::NRTDtorPP: 
			delete ((PolyPitchUGen*)cmd->nrtallocated);
			return true;
	}
	
	return false;
}

bool cmdStage3(World* world, CmdData* cmd) // RT
{
	switch (cmd->type) {
		case CmdData::NRTCtorPP: {
			((PolyPitch*)cmd->unit)->ugen_ = (PolyPitchUGen*)cmd->nrtallocated;
		}
			return true;
	}
	return false;
}

bool cmdStage4(World* world, CmdData* cmd) // NRT
{
	return true;
}

void cmdCleanup(World* world, void* cmd)
{
	RTFree(world, cmd);
}






//if just using for construction and destruction, can have named functions for construction and destruction, don't need to share same function with switch on type

void PolyPitch_Ctor( PolyPitch* unit ) {
	
    unit->maxvoices_ = ZIN0(1); 
    unit->levelcompressionfactor_= ZIN0(2); 
    unit->mixleftterm_= ZIN0(3);
    unit->torprec_= ZIN0(4);
    unit->cancellationweight_= ZIN0(5); 
    unit->polyphonyestimategamma_= ZIN0(6);

	unit->ugen_= 0; 
    
    //check sampling rate is 44100 and block size is 64
    unit->m_blocksize = unit->mWorld->mFullRate.mBufLength;
	 
	if(unit->m_blocksize!=64) { 
		printf("PolyPitch complains: block size not 64, you have %d\n", unit->m_blocksize);
		SETCALC(*ClearUnitOutputs);
		unit->mDone = true; 
		return; 
	}
	
	unit->m_sr = unit->mWorld->mSampleRate;
	
//	if(unit->m_sr!=44100) {
//        printf("PolyPitch complains: sample rate not 44100, you have %d\n", unit->m_sr);
//        SETCALC(*ClearUnitOutputs);
//        unit->mDone = true; 
//        return; 
//	}
//    

    //since common sample rate for all UGens running on this scsynth instance, shouldn't be any conflict
    g_samplingrate = unit->m_sr; 
    g_sampleperiod = 1.0f/g_samplingrate; 
    
    
	CmdData* cmd = (CmdData*)RTAlloc(unit->mWorld, sizeof(CmdData));
	cmd->samplingrate_ = unit->mRate->mSampleRate; 
	cmd->unit = (Unit *)unit; 
	cmd->nrtallocated = NULL; //will be allocated in NRT thread 
	cmd->type = CmdData::NRTCtorPP; 
	
	//(AsyncStageFn)PitchNoteUGencmdStage4
	
	DoAsynchronousCommand(unit->mWorld, 0, "", (void*)cmd,
						  (AsyncStageFn)cmdStage2,
						  (AsyncStageFn)cmdStage3,
						  NULL,
						  cmdCleanup,
						  0, 0);
	
		SETCALC(PolyPitch_next);
}





void PolyPitch_Dtor(PolyPitch *unit)
{
	
	if (unit->ugen_) {
		
		CmdData* cmd = (CmdData*)RTAlloc(unit->mWorld, sizeof(CmdData));
		cmd->unit = unit; 
		cmd->nrtallocated = unit->ugen_; 
		unit->ugen_ = NULL; //no longer available, will be deallocated in NRT thread
		cmd->type = CmdData::NRTDtorPP; 

		DoAsynchronousCommand(unit->mWorld, 0, "", (void*)cmd,
							  (AsyncStageFn)cmdStage2,
							  NULL,NULL,
							  cmdCleanup,
							  0, 0);
		
	}
	
}



//
//void PolyPitch_Ctor( PolyPitch* unit ) {
//	
//	float samplingrate = unit->mRate->mSampleRate; //unit->mWorld->mFullRate.mSampleRate; 
//	
//	SETCALC(PolyPitch_next);
//}
//
//void PolyPitch_Dtor(PolyPitch *unit)
//{
//
//	
//}




void PolyPitch_next( PolyPitch *unit, int inNumSamples ) {
	
    float *input = IN(0); 
    //float *output = OUT(0);
	
    int numSamples = unit->mWorld->mFullRate.mBufLength;    //blocksize should be 64
    
    //amortisation carried out internally? need to convert from ugen to current output format
    
    //unit->ugen_->compute( input, output, inNumSamples ); 
	
    PolyPitchUGen * pUGen = unit->ugen_; 
    
    if (pUGen != 0) {
        
        pUGen->compute(input,numSamples); 
        
        PeriodicityAnalysis * ppp = &pUGen->polypitch_;
        
        int voices = ppp->numvoicesdetected_; 
        
        ZOUT0(0) = voices; 
        
        //printf("numvoices %d ",voices); 
        
        //        cout << i << endl;
        for  (int j = 0; j<voices; ++j) {
            //            
            //            //cout << i << " " << polypitch.winningtor_ << " " << polypitch.bestsalience_ << endl; 
            //std::cout << j << " " << (1.0/polypitch_.voiceperiods_[j]) << " " << polypitch_.voicesaliences_[j] << std::endl; 
            //  
            
            int index= 2*j+1;
            
            ZOUT0(index) = 1.0f/(ppp->voiceperiods_[j]); 
            
            ++index; 
            
            ZOUT0(index) = ppp->voicesaliences_[j];
            
        }
        
        
    }
    
    else {
		
		ZOUT0(0) = 0.0; //no voices detected
	}
    
    
//    for (int j=0; j<inNumSamples; ++j)
//        out[j]= 0.0; 
    
}









PluginLoad(PolyPitch) {
	
	ft = inTable;
	
    DefineDtorCantAliasUnit(PolyPitch);
	
	
}





