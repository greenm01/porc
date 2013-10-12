#PORC Setup Guide

A quick tutorial on how to get up and running for PORC and other Python scientific 
tools for speaker and room design/modeling tools

##Configure Python - All OS


anaconda python distribution - https://store.continuum.io/cshop/anaconda/

and instructions on http://docs.continuum.io/anaconda/install.html for the install process.

##Install libsndfile - Windows



##Install libsndfile - Linux

On Linux install libsndfile-dev from your Linux distribution package management (apt-get, 
yum etc)

##Install libsndfile - OS X

Use Homebrew http://brew.sh/ as package management for OS X, many more tools available in 
this almost indispensable utility. Follow install instructions on Homebrew web page. 
One-line install process at terminal prompt.

Once Homebrew is setup run:

    brew install libsndfile
    
##Install scikits.audiolab

With Python and dependencies set up you should now have pip working so:
    pip install scikits.audiolab

##Download and run PORC
Download PORC - https://github.com/zzzzrrr/porc/archive/master.zip

Unzip the files to a new folder. Now at terminal prompt navigate to the directory where 
PORC files are stored.

    python porc.py -h 

Now just follow instructions in the help file.