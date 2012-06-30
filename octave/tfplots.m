
function tfplots(data, color, Fs, fract);

%TFPLOT - Smoothed transfer fucntion plotting
%   TFPLOTS(IMPRESP,COLOR, Fs, FRACT)
%   Logarithmic transfer function plot from impluse response IMPRESP. 
%   A half hanning window is applied before a 2^18 point FFT, then the data is colleced
%   into logaritmically spaced bins and the average power is computed for
%   each bin (100/octave). Then this is power-smoothed by a hanning window, where
%   FRACT defines the fractional-octave smoothing (default is 3, meaning third-octave).
%   The length of the smoothing hanning window is the double compared to the distance
%   defined by FRACT.
%   The sampling frequency is set by FS (default is 44.1 kHz) and the plotting color is set by the COLOR variable
%   (default is 'b').
%
%   C. Balazs Bank, 2007.

octbin=100;

if nargin==1,
    color='b';
    fract=3;
    Fs=44100;
end;

if nargin==2,
    fract=3;
    Fs=44100;
end;

if nargin==3,
    fract=3;
end;

FFTSIZE=2^18;
data=data(:);

logfact=2^(1/octbin);
LOGN=floor(log(Fs/2)/log(logfact));
logscale=logfact.^[0:LOGN]; %logarithmic scale from 1 Hz to Fs/2

WL=length(data);
hann=hanning(WL*2);
endwin=hann(WL+1:2*WL);
tf=fft(data.*endwin,FFTSIZE);

magn=(abs(tf(1:FFTSIZE/2)));
compamp=tf(1:FFTSIZE/2);


%creating 100th octave resolution log. spaced data from the lin. spaced FFT data
clear logmagn;
fstep=Fs/FFTSIZE;
for k=0:LOGN,
	start=round(logscale(k+1)/sqrt(logfact)/fstep);
	start=max(start,1);
	start=min(start,FFTSIZE/2);
	stop=round(logscale(k+1)*sqrt(logfact)/fstep);
	stop=max(stop,1);
	stop=min(stop,FFTSIZE/2);
	logmagn(k+1)=sqrt(mean(magn(start:stop).^2)); %averaging the power
end;

%creating hanning window
HL=2*round(octbin/fract); %fractional octave smoothing
hh=hanning(HL);

L=length(logmagn);
logmagn(L+1:L+HL)=0;

%Smoothing the log. spaced data by convonvling with the hanning window
tmp=fftfilt(hh,logmagn.^2);
smoothmagn=sqrt(tmp(HL/2+1:HL/2+L)/sum(hh));

%plotting
semilogx(logscale,20*log10(smoothmagn),color);
