/*
 *
 *	csound-xtract demonstration
 *	By Richard Knight 2021
 *
 *	This demonstrates some of the csound-xtract plugin opcodes in three instruments:
 *	
 *	1. test_dump		prints selected analysed features of a sound
 *	2. test_match		uses feature matching to play a corpus based on analysis of an input sound
 *	3. test_distance	prints the distance between two streams of analysed features
 *
 *	
 *	csound-xtract is current experimental. Using different buffer and block sizes for xtprofile may
 *	result in quite different feature results depending on the input sound. For certain opcodes, namely
 *	xtcorpusmatch, this require some fine-tuning to obtain the desired results based on the character of 
 *	the input sound.
 *
 */

<CsoundSynthesizer>
<CsOptions>
-odac
</CsOptions>
<CsInstruments>
sr = 44100
kr = 4410
nchnls = 2
0dbfs = 1
seed 0


; source and corpus sounds
gSsource = "sounds/fox.wav"
gScorpus = "sounds/magicshop.wav"

; corpus f-table
gicorpuswave ftgen 0, 0, 0, 1, gScorpus, 0, 0, 1

; half size window for sndwarp reading
giwindow ftgen 0, 0, 512, 9, 0.5, 1, 0 




/*
	play the source sound, dry, and print the extracted features
*/
instr test_dump
	; use centroid, irregularity, sharpness and smoothness
	iprofile xtprofile 2048, 1024, 0, 1, 0, 0, 0, 1, 0, 1, 1

	; read the source/driving sound
	asource diskin gSsource

	; perform feature extraction on source sound
	ixtract, kdone xtractor iprofile, asource

	; print the features
	if (kdone == 1) then
		kanalysis[] xtdump ixtract
		printarray kanalysis
	endif

	outs asource, asource
endin



/* 
	match the source sound against the corpus and play the 'best match'
*/
instr test_match

	; use MFCC, rms, flatness, irregularity, spectral power, sharpness and smoothness
	iprofile xtprofile 2048, 1024, 1, 0, 0, 1, 1, 1, 1, 1, 1

	; analyse the corpus wave
	icorpus xtcorpus iprofile, gicorpuswave

	; read the source/driving sound
	asource diskin gSsource
	
	; perform feature extraction on source sound and match with corpus
	ixtract, kdone xtractor iprofile, asource
	kdex xtcorpusmatch icorpus, ixtract, kdone

	; if matched index != 0 then calculate the position to be used for sndwarp
	if (kdex != 0) then
		ilen = ftlen(gicorpuswave)
		icsr = ftsr(gicorpuswave)
		icduration = ilen / icsr
		icps = 1/(icduration)
		aphs, a_ syncphasor icps, a(kdone)
		apos = (((aphs * ilen) + kdex) / ilen) * icduration

	endif

	; play the position in the corpus wave
	amatched sndwarp 0.7, apos, 1, gicorpuswave, 0, 1024, 512, 8, giwindow, 1	
	outs amatched, amatched
endin



/*
	print the distance between two sound streams
*/
instr test_distance
	; use all features
	iprofile xtprofile 2048, 1024

	; read the source/driving sound
	asource diskin gSsource

	; modify
	amodified = butterlp(asource, 10000)

	; perform feature extraction on sounds
	ixtract1, kdone1 xtractor iprofile, asource
	ixtract2, kdone2 xtractor iprofile, amodified

	; calculate the distance and print it
	kdistance xtdistance ixtract1, ixtract2, 1
	printk2 kdistance
endin



</CsInstruments>
<CsScore>

i"test_dump" 0 3
i"test_match" 3 3
i"test_distance" 6 3

</CsScore>
</CsoundSynthesizer>