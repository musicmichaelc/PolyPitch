//live analysis of frequency track data

b= Buffer.read(s,"/data/audio/mirdata/scmirtests/pianotest.wav");

b= Buffer.read(s,"/data/audio/mirdata/scmirtests/pixiesivebeentired.wav");

b= Buffer.read(s,"/data/audio/mirdata/scmirtests/Everybody Loves Somebody.wav");



(
SynthDef(\pptest1,{arg amp=0.1, freq=440, dur=0.1;

var env, sound;

env = EnvGen.ar(Env([0,1,1,0],[0.01,dur,0.01]),doneAction:2);

sound = LFCub.ar(freq+(LFNoise1.kr(Rand(5,7.5),Line.kr(0.0,Rand(2.0,5.0),dur))));

Out.ar(0,Pan2.ar(sound*env*amp,0.0));

}).add;

SynthDef(\pptest2,{arg amp=0.1, freq=440, dur=0.1;

var env, sound;

env = EnvGen.ar(Env([0,1,1,0],[0.01,dur,0.01]),doneAction:2);

sound = LFSaw.ar(freq+(LFNoise1.kr(Rand(5,7.5),Line.kr(0.0,Rand(2.0,5.0),dur))));

Out.ar(0,Pan2.ar(sound*env*amp,0.5));

}).add
)

Synth(\pptest1,[\amp, 0.1, \freq, 440, \dur, 0.5])


//better for piano polyphony
(
{

var in, tracked;

//in= SoundIn.ar;
in = PlayBuf.ar(1,b,loop:1);

tracked=PolyPitch.kr(in,5,-0.1,4.0,cancellationweight:0.5); //0.0000001,0.5,0.9

SendReply.kr(Impulse.kr(44100/2048),'polyparam',tracked);

Out.ar(0,in);

}.play
)


//better for complex polyphonic audio with background drums?
(
{

var in, tracked;

//in= SoundIn.ar;
in = PlayBuf.ar(1,b,loop:1);

tracked=PolyPitch.kr(in,10,-0.1,10.0,cancellationweight:1.0); //0.0000001,0.5,0.9

SendReply.kr(Impulse.kr(44100/2048),'polyparam',tracked);

Out.ar(0,in);

}.play
)




(
var maxvoices = 10;
var activevoices = List[];
var prevnotes = 0.0!maxvoices;
var threshold = 0.3; //1/5 of a semitone resolution to say frequency is the 'same'. Vibrato can still cause a confound, though unlikely to move too far (some singers might show 1 ST deviation at 7 Hz for example, versus 20Hz sampling rate of PolyPitch, though vibrato is smeared by the FFTs)
var threshold2 = 0.4; //account for wider vibrato over longer term

o.remove;

o = OSCresponder(nil, 'polyparam', { |t, r, msg|

	var numvoices = msg[3];
	var newnotes = List[];

	numvoices.do{|i|
		var indexnow = 4+(i*2);
		var note = msg[indexnow].cpsmidi;
		var salience = msg[indexnow+1]; //will be used for strength of trails

		//[note,salience].postln;

		Synth(\pptest1,[\amp, (0.3*salience).min(0.1), \freq, 2.0*msg[indexnow], \dur, 0.1]);

		block {|break|

		prevnotes.do{|notethen|

			//octave errors mess this up, can hear some in general response that are picking up pitch content that
			//isn't octave consistent frame to frame

			var diff = abs(note-notethen);
			//var voicealready= false;

			if(diff<threshold) {

				//create active voice if not one close already; else add to salience of active voice match

				activevoices.do{|voice|
					var voicenote = voice[0];
					var voicediff = abs(note-voicenote);

					if(voicediff<threshold2) {

						voice[1] = voice[1] + salience;

						break.value;
					}

				};

				//if got to here, no existing matching voice

				//third param is the survival frames counter
				activevoices.add([note, salience,0] );

				break.value;
			};

			};

		};

		newnotes.add(note);
	};

	//update active voices, kill off the weak, reduce all ratings
	activevoices.do {|voice|   voice[1] = voice[1] * 0.8; voice[2] = voice[2]+1;
		if(voice[2]==1,{Synth(\pptest2,[\amp, 0.02, \freq, voice[0].midicps, \dur, 0.2+(0.5*(voice[1].min(1.0)))])})
	};

	activevoices = activevoices.select {|voice|   voice[1] > 0.25 };

	activevoices.postln;

	prevnotes = newnotes;

}).add;

)


