//live analysis of frequency track data

(
{

var in, tracked;

in= SoundIn.ar;

tracked=PolyPitch.kr(in,10);

SendReply.kr(Impulse.kr(44100/2048),'polyparam',tracked);

}.play
)



//12TET chromatic
(
var maxvoices = 5;
var chromascores = 0.0!12;
var leak = 0.9;
var ms;

w=Window.new("PolyPitch Chromagram", Rect(200,400,300,300));

ms= MultiSliderView.new(w, Rect(10,10,240,280));

ms.value_(Array.fill(12,0.0));
ms.valueThumbSize_(20.0);
ms.indexThumbSize_(20.0);
ms.gap_(0);

w.front;

r.stop;

o.remove;

o = OSCresponder(nil, 'polyparam', { |t, r, msg|

	var numvoices = msg[3];
	var newnotes = List[];

	numvoices.do{|i|
		var indexnow = 4+(i*2);
		var chroma = (round(msg[indexnow].cpsmidi,1.0))%12;
		var salience = msg[indexnow+1]; //will be used for strength of trails

		chromascores[chroma] = chromascores[chroma] + salience;


	};

	chromascores = chromascores*leak;

}).add;



r= {

inf.do{

{ms.value_(chromascores*0.2)}.defer;

(2048.0/s.sampleRate).wait; //25 frames per second
};

}.fork;


)





//MIDI Note gram (24 to 95)
(
var maxvoices = 5;
var midiscores = 0.0!72;
var leak = 0.9;
var ms;

w=Window.new("PolyPitch MIDI Note gram (MIDI 24 to 95)", Rect(200,400,300,300));

ms= MultiSliderView.new(w, Rect(10,10,288,280));

ms.value_(Array.fill(72,0.0));
ms.valueThumbSize_(4.0);
ms.indexThumbSize_(4.0);
ms.gap_(0);

w.front;

r.stop;

o.remove;

o = OSCresponder(nil, 'polyparam', { |t, r, msg|

	var numvoices = msg[3];
	var newnotes = List[];

	numvoices.do{|i|
		var indexnow = 4+(i*2);
		var notenow = round(msg[indexnow].cpsmidi,1.0);
		var chroma = (notenow-24)%72;
		var salience = msg[indexnow+1]; //will be used for strength of trails

		if ((notenow<24) || (notenow>=96)) {salience = 0.0};

		midiscores[chroma] = midiscores[chroma] + salience;


	};

	midiscores = midiscores*leak;

}).add;



r= {

inf.do{

{ms.value_(midiscores*0.2)}.defer;

(2048.0/s.sampleRate).wait; //25 frames per second
};

}.fork;


)



