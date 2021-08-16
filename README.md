# csound-xtract : Csound feature extraction using libXtract

## Overview
csound-xtract is a set of plugin opcodes which use libXtract to perform feature extraction and associated tasks from within Csound.

Development is still ongoing and subject to various research matters, thus is provided in an experimental/alpha state and may contain bugs. Parts of the code are due overhauls and refactoring, but the intention is for the opcodes and general operation to remain the same as presented here.


## Requirements
* Cmake >= 2.8.12
* Csound with development headers >= 6.14.0
* [LibXtract](https://github.com/jamiebullock/LibXtract)

Tested on Linux and Windows 7 with MSYS as of March 2021.


## Installation
Create a build directory at the top of the source tree, execute *cmake ..*, *make* and optionally *make install* as root. If the latter is not used/possible then the resulting libcsxtract library can be used with the *--opcode-lib* flag in Csound.
eg:

	mkdir build && cd build
	cmake ..
	make && sudo make install

Cmake should find Csound and libXtract using the modules in the cmake/Modules directory and installation should be as simple as above.

## Examples
Some examples are provided in the examples directory.


## Opcode reference

### iprofile xtprofile [ibuffersize=4096, iblocksize=512, imfccs=1, icentroid=1, izerocrossings=1, irms=1, iflatness=1, iirregularity=1, ipower=1, isharpness=1, ismoothness=1]
Declare a feature extraction profile for use in extraction opcodes. The exact techniques for extraction of individual features can be found by examining the libXtract documentation and source code.

* iprofile : the profile handle

* ibuffersize : buffer size used in extraction
* iblocksize : block size for extraction
* imfccs : use MFCCs
* icentroid : use spectral centroid
* izerocrossings : use zero crossings
* irms : use RMS
* iflatness : use spectral flatness
* iirregularity : use spectral irregularity
* ipower : use spectral power
* isharpness : use spectral sharpnesss


### icorpus xtcorpus iprofile, ifn
Analyse sound contained in a f-table and store in a handle, for later use in comparison/matching opcodes. Done during init time.

* icorpus : the corpus handle

* iprofile : profile handle as created by xtprofile
* ifn : f-table containing the sound to be analysed, typically GEN1


### ixtract, kdone xtractor iprofile, ain
Analyse a live sound.

* ixtract : the extraction handle

* iprofile : profile handle as created by xtprofile
* ain : sound to analyse


### kdistance xtdistance ixtract1, ixtract2, ktrigger, [idistancefunc=0]
Compare two extraction streams using a basic distance function between each frame of the analyses. The profiles used for the streams must be the same.

* kdistance : the calculated distance, 0 should represent no difference between the analyses.

* ixtract1 : analysis stream as created by xtractor
* ixtract2 : analysis stream as created by xtractor
* ktrigger : comparison is conducted when 1
* idistancefunc : 0 for Euclidean distance, 1 for Manhattan distance


### kanalysis[] xtdump ixtract
Obtain the analysis data from an xtractor handle. The array length will depend on the number of features specified in the profile used (MFCC uses 13 indexes, all others use 1). The indexes are presented in the same order as declared in xtprofile. For example, a profile using MFCC and centroid would mean indexes 0 to 12 would be MFCCs, and 13 would be centroid. Similarly a profile using only centroid and flatness would imply index 0 would be centroid, and 1 would be flatness.

* kanalysis[] : the analysed features

* ixtract : analysis stream as created by xtractor


### kdone, kanalysis[] xtaccdump ixtract, ktrigger
Obtain the analysis data from an xtractor handle as with xtdump, but accumulate the analyses and output the mean when ktrigger is 1.

* kdone : outputs 1 when new data is provided, 0 at all other times
* kanalysis[] : the analysed features as with xtdump. 

* ixtract : analysis stream as created by xtractor
* ktrigger : triggers the output of data when 1


### kposition xtcorpusmatch icorpus, ixtract, ktrigger, [idistancefunc=0]
Obtain the nearest match in a corpus created with xtcorpus by comparing a live input stream created by xtractor.

* kposition : position of nearest corpus region, in sample points

* icorpus : corpus handle as created by xtcorpus
* ixtract : analysis stream as created by xtractor
* ktrigger : perform the comparison when 1
* idistancefunc : 0 for Euclidean distance, 1 for Manhattan distance



