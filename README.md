[![Build Status](https://travis-ci.com/randomsamples/porc.svg?branch=master)](https://travis-ci.com/randomsamples/porc)

Python Open Room Correction (PORC)
==================================

PORC now includes mixed-phase compensation (see below)!

DSP Loudspeaker-Room correction filter wizard; transfer function modeling and equalization 
by fixed-pole parallel filters. Algorithm ported to Python by Mason A. Green, based on the work 
of Dr. Balazs Bank: http://home.mit.bme.hu/~bank/parfilt/

More details about Dr. Bank's parallel filter can be found in the papers:

    Balazs Bank, "Perceptually Motivated Audio Equalization Using Fixed-Pole Parallel
    Second-Order Filters," IEEE Signal Processing Letters, 2008.

http://www.acoustics.hut.fi/go/spl08-parfilt

    Balazs Bank, "Direct Design of Parallel Second-order Filters for
    Instrument Body Modeling," International Computer Music Conference,
    Copenhagen, Denmark, Aug. 2007.

http://www.acoustics.hut.fi/go/icmc07-parfilt
    
Mixed-Phase Compensation references:

	Alberto Carini, et al, "Mixed Time-Frequency approach for Multipoint
	Room Response Equalization," AES 45th International Conference, 2012
	
	Defrance & Polak, "Measuring the mixing time in auditoria," Acoustics
	Paris 2008
	
Required Python dependencies:

    1) Python 2.7
    2) Scientific Python: SciPy, Numpy, & Matplotlib

The easiest install method on Windows is simply to install the continuum.io Anaconda package.

Measurement
===========

One needs to measure the log-frequency impulse response of your speakers with a 
calibrated Electret Measurement Microphone, e.g. Dayton Audio EMM-6. Software 
such as Room EQ Wizard (REQ), Holm Impulse, or Arta may be used for this purpose:
http://www.hometheatershack.com/roomeq/

Usage
=====

porc.py [-h] [--mixed] [-t FILE] [-n NTAPS] [-o OPFORMAT] input_file output_file

    python porc.py -t tact30f.txt -n 6144 -o bin l48.wav leq48.bin

Use the -h flag for help!

PORC has been tested successfully on both Linux and Windows 7 with Python 2.7. Linux depenency 
install is fairly straightforward. Windows install packages are available for all dependencies.

Target Response
===============

The default target curve for PORC is flat. Included in the data directory are a number 
of target curves. Experiment to suit your listening preferences. Use the [-t] flag to load a target
file.

One may also target a flat curve, and then use separate parametric equalization for bass boosting
and other pschoaccoustic preferences. 

For further reference, the B&K House Curve is a good place to start. Read "Relevant loudspeaker 
tests in studios in Hi-Fi dealers' demo rooms in the home etc.," Figure 5:
http://www.bksv.com/doc/17-197.pdf

Mixed-Phase Compensation
==============

To use mixed-phase compensation, one needs to specify the [--mixed] flag. One also needs to modify
the Room Impulse Response (RIR) to remove leading silence (zeros) before the main impulse. You can
easily do this with Audacity or REQ.

Example:

	python porc.py --mixed -t tact30f.txt -n 6144 -o bin l48.wav leq48.wav
	
Have some patience with this method. The convolution takes a few CPU cycles.

PC Convolution
==============

Suggestions:

Windows (foobar2000 convolver)
Linux (jconvolver w/ jcgui & Jack)

You may need to merge left and right channels into a single stereo .wav 

    sox -M le148.wav req48.wav equalizer.wav


OpenDRC Convolution
===================

Use -o bin flag to set binary 32bit IEEE floating point mono file format output for OpenDRC.


TODO
====

	Implement algo to automatically remove leading silence (zeros) from RIR.
	Add a GUI Frontend (pretty interactive graphs, drawing target curve, etc...)
	Update this page with better documentation!


Contact
=======

Complaints, suggestions, bugfixes: mason dot green at gmail
