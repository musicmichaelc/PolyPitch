#PolyPitch
##A [SuperCollider](https://supercollider.github.io/) plugin
This is a SC plug-in for multiple fundamental frequency tracking, after Anssi Klapuri\'s 2008 paper \'Multipitch analysis of polyphonic music and speech signals using an auditory model\'


PolyPitch is released under [GNU GPL 3](sources/COPYING), for [SuperCollider 3](https://supercollider.github.io/), by [Nick Collins](https://composerprogrammer.com/index.html)


###Installation
To install: put the PolyPitch folder in your SC extensions directory. Copy into that folder a PolyPitch.scx file. There are some prebuilt ones provided, and compilation instructions also follow. 

Prebuilt binaries can be found on the "release" tab above.

###Compiling
To compile for platforms other than Mac, using the SC native FFT interface, comment out:
#define OSXACCELERATENOTSCFFT
in MultipitchHeader.h

Compilation requires `libsamplerate` from http://www.mega-nerd.com/SRC/

####Compilation example: 64 bit \(e.g. SC 3.6\)
Make sure in `CMakeLists.txt` that the 64bit library linking line is uncommented and 32bit is commented
Then on the command line \(assuming `CMake` is installed and you have the source code of SuperCollider at `/data/gitprojects/scdev/supercollider`\):
```
cmake -DSC_PATH=/data/gitprojects/scdev/supercollider -DCMAKE_OSX_ARCHITECTURES='x86_64'
make
```

####Compilation example: 
32 bit on OS X

Make sure in `CMakeLists.txt` that the 32bit library linking line is uncommented and 64bit is commented

Then on the command line \(assuming `CMake` is installed and you have the source code of SuperCollider at /data/gitprojects/supercollider\):
```
cd into build directory
cmake -DSC_PATH=/data/gitprojects/supercollider -DCMAKE_OSX_ARCHITECTURES='i386' ..
make
```
 
