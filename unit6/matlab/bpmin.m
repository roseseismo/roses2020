function epb=bpmin(H,F,dB)
% epb=BPMIN(H,F,dB)
%
% INPUT:
%
% H     Magnitude response (in deciBel)
% F     Frequency vector
% dB    Stopband level (positive)
%
% OUTPUT:
% 
% epb   Frequencies where the pass band exceeds dB level
%
% Last modified by fjsimons-at-alum.mit.edu, July 23rd, 2003

defval('dB',3)

H=H+dB;
Hm=find(H==max(H));

% Both are now monotonic functions
F1=interp1(H(1:Hm),F(1:Hm),0);
F2=interp1(H(Hm+1:end),F(Hm+1:end),0);

epb=[F1 F2];


