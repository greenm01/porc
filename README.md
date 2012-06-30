Python Open Room Correction (PORC)
==================================

DSP Loudspeaker-Room correction filter wizard; transfer function modeling and equalization 
by fixed-pole parallel filters.

Algorithm ported to Python by Mason A. Green, based on the work of Dr. Balazs Bank:
http://home.mit.bme.hu/~bank/parfilt/

More details about Dr. Bank's parallel filter can be found in the papers:

	Balazs Bank, "Perceptually Motivated Audio Equalization Using Fixed-Pole Parallel
	Second-Order Filters", IEEE Signal Processing Letters, 2008.
	http://www.acoustics.hut.fi/go/spl08-parfilt

	Balazs Bank, "Direct Design of Parallel Second-order Filters for
	Instrument Body Modeling", International Computer Music Conference,
	Copenhagen, Denmark, Aug. 2007.
	http://www.acoustics.hut.fi/go/icmc07-parfilt
	
To run PORC, you need to install the following Python libs:

	1) SciPy
	2) Numpy
	3) Matplotlib
	
One needs to measure the log-frequency impulse response of your speakers with a 
calibrated Electret Measurement Microphone such as the Dayton Audio EMM-6. Software 
such as Room EQ Wizard (REQ) may be used for this purpose:
http://www.hometheatershack.com/roomeq/

TODO:

	Update this page with detailed instructions
	Explain OpenDRC integration
	Add target splines

