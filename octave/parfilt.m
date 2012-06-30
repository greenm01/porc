
function y=parfilt(Bm,Am,FIR,x);

%PARFILT filtering by parallel second-order sections.
%   Y=PARFILT(Bm,Am,FIRcoeff,X) filters the signal X by the parallel
%   second-order sections given in the matrices Bm and Am. Additionally, X
%   is filtered by an FIR filter given by FIRcoeff, and this is also added
%   to the output. If there is no parallel FIR part (not even a simple
%   coefficient), then use FIRcoeff=0 as a parameter.
%
%   The Bm and Am matrices are containing the [b0 b1]' and [1 a0 a1]'
%   coefficients for the different sections in their columns. For example,
%   Bm(:,3) gives the [b0 b1]' parameters of the third second-order
%   section. The parallel second-order sections can be designed by the
%   PARFILTDES command.
%
%   More details about the parallel filter can be found in the papers
%
%	 Balazs Bank, "Perceptually Motivated Audio Equalization Using Fixed-Pole Parallel
%   Second-Order Filters", IEEE Signal Processing Letters, 2008.
%   http://www.acoustics.hut.fi/go/spl08-parfilt
%
%   Balazs Bank, "Direct Design of Parallel Second-order Filters for
%   Instrument Body Modeling", International Computer Music Conference,
%   Copenhagen, Denmark, Aug. 2007.
%   http://www.acoustics.hut.fi/go/icmc07-parfilt
%
%   C. Balazs Bank, Helsinki University of Technology, 2007.

y=zeros(size(x));
s=size(Am);
for k=1:s(2),
    y=y+filter(Bm(:,k),Am(:,k),x);
end;
y=y+filter(FIR,1,x);
    

