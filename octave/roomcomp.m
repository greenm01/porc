#! octave 
% Adapted from Balazs Bank's MATLAB code for OpenDRC
% Modified by Mason A. Green June 2012
%
% Balazs Bank: Perceptually Motivated Audio Equalization Using Fixed-Pole Parallel Second-Order Filters
% IEEE Signal Processing Letters, 2008
%
% C. Balazs Bank, Helsinki University of Technology, 2007.

arg_list = argv ();

if length(arg_list) != 2;
   printf("\nIncorrect command line arguments!\n");
   printf ("\n Format: octave roomcomp.m <impulse response input> <filter output>\n");
   printf("Example: octave roomcomp.m left_48khz.wav left_filter.wav\n");
   exit();
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logarithmic pole positioning

[data, Fs] = auload(arg_list{1});

% You may need to change R depending on the number of poles 
% (more poles: larger, less poles: smaller)
R = 0.5; 
% Two sets of log. resolution
fplog=[logspace(log10(30),log10(200),13) logspace(log10(250),log10(18000),12)]; 

wp=2*pi*fplog/Fs;
p=R.^(wp/pi).*exp(j*wp); 
plog=[p conj(p)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing data

[cp,minresp]=rceps(data); %making the measured response minumum-phase
output=zeros(length(minresp),1);
output(1)=1; %target 

[Bf,Af]=butter(4,30/(Fs/2),'high');
outf=filter(Bf,Af,output); %making the target output a 30 Hz highpass

imp=zeros(1,length(data));
imp(1)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter design

%Parallel filter design
[Bm,Am,FIR]=parfiltid(minresp,outf,plog,1); 

% equalized loudspeaker response - filtering the 
% measured transfer function by the parallel filter
equalizedresp=parfilt(Bm,Am,FIR,data); 
									   
% equalizer impulse response - filtering a unit pulse
equalizer=parfilt(Bm,Am,FIR,imp); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize: output filter and plot

% Convert to the time domain
y = ifft(equalizer);
% Downsample by n
n = ceil(length(y)/6148);
y = decimate(y, n);
% Save output filter
yy = fft(y);
yy = yy/max(abs(yy));
ausave(arg_list{2}, yy, Fs);

% Plot data. 
% NOTE: You may need to adjust the scaling factors to make 
%       your data fit nicely on the graph

%original loudspeaker-room response
tfplot(50*data,'b', Fs, 100, 'abs'); 
hold on; 
tfplots(50*data, 'r--', Fs); %3rd octave smoothed

%equalized loudspeaker-room response
tfplot(equalizedresp*0.12,'b', Fs, 100, 'abs'); 
tfplots(equalizedresp*0.12,'r--', Fs); %3rd octave smoothed

%equalizer transfer function
tfplot(yy,'b', Fs); 
L=[-2;2]*ones(1,length(fplog));	%indicating pole frequencies
line([fplog;fplog],L,'color','k');
hold off;

% Label plot
axis([20 20000 -60 60]);
set(gca,'FontName','Times','Fontsize',14);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');

text(500,25,'Unequalized loudspeaker-room response','FontName','Times','FontSize',10);
text(1000,14,'Equalizer transfer function','FontName','Times','FontSize',10);
text(1000,10,'(Black lines: pole locations)','FontName','Times','FontSize',10);
text(200,-43,'Equalized loudspeaker-room response','FontName','Times','FontSize',10);
print("eqplot.png","-dpng");

% Output messages
printf("\nGraphical plots written to eqplot.png\n"); 
printf("\nOutput filter length = %i taps\n", length(yy));
printf("Output filter written to %s\n", arg_list{2});

printf("\nUse sox to convert output wav to raw 32 bit IEEE floating point\n");
printf("Example: sox filter.wav -t f32 filter.bin\n");
