Python Open Room Correction (PORC)
==================================

DSP Loudspeaker-Room correction filter wizard; transfer function modeling and equalization 
by fixed-pole parallel filters. Algorithm ported to Python by Mason A. Green, based on the work 
of Dr. Balazs Bank: http://home.mit.bme.hu/~bank/parfilt/

More details about Dr. Bank's parallel filter can be found in the papers:

	Balazs Bank, "Perceptually Motivated Audio Equalization Using Fixed-Pole Parallel
	Second-Order Filters", IEEE Signal Processing Letters, 2008.

http://www.acoustics.hut.fi/go/spl08-parfilt

	Balazs Bank, "Direct Design of Parallel Second-order Filters for
	Instrument Body Modeling", International Computer Music Conference,
	Copenhagen, Denmark, Aug. 2007.

http://www.acoustics.hut.fi/go/icmc07-parfilt
	
Required Python dependencies:

	1) Scientific Python: SciPy, Numpy, & Matplotlib
	2) Scikits.audiolab
	3) libsndfile

This is a command line tool. Matplotlib will produce very nice graphs; if you desire a GUI 
feel free to fork this project.

Measurement
===========

One needs to measure the log-frequency impulse response of your speakers with a 
calibrated Electret Measurement Microphone, e.g. Dayton Audio EMM-6. Software 
such as Room EQ Wizard (REQ) may be used for this purpose:
http://www.hometheatershack.com/roomeq/

Usage
=====

porc.py [-h] [-t FILE] [-n NTAPS] I F

	python porc.py -t tact30f.txt -n 6148 l48.wav leq48.wav

Use the -h flag for help!

PORC has been tested successfully on both Linux and Windows 7 with Python 2.7. Linux depenency 
install is fairly straightforward. Windows install packages are available for all dependencies.

Target Response
===============

The default target curve for PORC is flat. Included in the data directory are a number 
of target curves. Experiment to suit your listening preferences (I prefer tact30f.txt, 
bk-48.txt).

For further reference, the B&K House Curve is a good place to start. Read "Relevant loudspeaker 
tests in studios in Hi-Fi dealers' demo rooms in the home etc.," Figure 5:
http://www.bksv.com/doc/17-197.pdf

PC Convolution
==============

Suggestions:

Windows (foobar2000 convolver)
Linux (jconvolver w/ jcgui & Jack)

You may need to merge left and right channels into a single stereo .wav 

    sox -M le148.wav req48.wav equalizer.wav


OpenDRC Convolution
===================

Use sox to convert output .wav to raw 32 bit IEEE floating point mono for the left & right channels

	sox leq48.wav -t f32 leq48.bin
	sox req48.wav -t f32 req48.bin

TODO
====

	Update this page with better documentation!
	
