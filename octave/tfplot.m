
function tfplot(data, color, Fs, octbin, avg);

%TFPLOT - Logarithmically smoothed transfer fucntion plot
%   TFPLOT(IMPRESP, COLOR, FS, OCTBIN, AVG)
%   Logarithmic transfer function plot from impluse response IMPRESP. 
%   A half hanning window is applied before a 2^18 point FFT, then the 
%   data is colleced into logaritmically spaced bins and the average 
%   response is computed for each bin. OCTBIN sets the number of bins 
%   in one octave, the default is 100 (lower numbers mean more smoothing). 
%   The sampling frequency is set by FS (default is 44.1 kHz) and the 
%   plotting color is set by the COLOR variable (default is 'b').
%
%   If the AVG variable is set to 'power' then the power is averaged
%   in the logaritmic bin, if it is 'abs' then the absolute value. If the
%   AVG parameter is set to 'comp' or omitted, it averages the complex
%   magnitude (i.e., this is the default).
%
%   C. Balazs Bank, 2006-2007.


if nargin<5,
    avg='comp';
end;

if nargin<4,
    octbin=100;
end;

if nargin<3,
    Fs=44100;
end;

if nargin<2,
    color='b';
end;


FFTSIZE=2^18;
data=data(:);

logfact=2^(1/octbin);
LOGN=floor(log(Fs/2)/log(logfact));
logscale=logfact.^[0:LOGN]; %logarithmic scale from 1 Hz to Fs/2

WL=length(data); %creating a half hanning window
hann=hanning(WL*2);
endwin=hann(WL+1:2*WL);
tf=fft(data.*endwin,FFTSIZE); %FFT
compamp=tf(1:FFTSIZE/2);

fstep=Fs/FFTSIZE;
for k=0:LOGN,
   %finding the start and end positions of the logaritmic bin
   
   start=round(logscale(k+1)/sqrt(logfact)/fstep); 
   start=max(start,1);
   start=min(start,FFTSIZE/2);
   stop=round(logscale(k+1)*sqrt(logfact)/fstep)-1;
   stop=max(stop,start);
   stop=max(stop,1);
   stop=min(stop,FFTSIZE/2);
   
   if strcmpi(avg,'comp'),   %averaging the complex transfer function
       logmagn(k+1)=abs(mean(compamp(start:stop))); 
   end;
   if strcmpi(avg,'abs'),
      logmagn(k+1)=mean(abs(compamp(start:stop))); 
   end;
   if strcmpi(avg,'power'),
       logmagn(k+1)=sqrt(mean(abs(compamp(start:stop)).^2)); 
   end;

   
end;


%figure;
semilogx(logscale,20*log10(logmagn),color);
