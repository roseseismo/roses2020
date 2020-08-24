function [xf,co,npol,npas,tipe,HABS2,F,EPB]=...
    bpass(x,Fs,colo,cohi,npol,npas,tipe)
% [xf,HABS2,F,EPB]=BANDPASS(x,Fs,colo,cohi,npol,npas,tipe)
%
% Filters signal 'x' with filter 'tipe' and corner
% frequencies 'cohi' and 'cohi' in Hz with 'npol' the 
% number of poles and in 'npas' passes. Sampling rate is 'Fs' in Hz.
%
% INPUT:
%
% x         The signal
% Fs        Its sampling frequency
% colo      The lower corner frequency
% cohi      The higher corner frequency
% npol      The number of poles
% npas      The number of passes
% tipe      The filter name
%
% OUTPUT:
%
% xf      The filtered signal
%
% Compare in SAC bp butter co 0.05 5 n 2 p 1
%
% Returns the npas frequency response and the effective pass band for
% one or two passes (3 dB level)
%
% You'll see that plot(F,decibel(HABS2)) (this is what freqz plots)
% shows how' you concentrate between cohi and colo at the 
% 3 dB-level
%
% Last modified by fjsimons-at-alum.mit.edu, Nov 21th, 2004

defval('npol',2)
defval('npas',1)
defval('colo',0.05)
defval('cohi',0.50)
defval('Fs',110)
defval('tipe','butter')

%disp(sprintf('BANDPASS %3.3f-%3.3f Hz %i pass %i poles %s',...
%	     colo,cohi,npas,npol,tipe))
						
% Corner frequency is in Hertz, now it is as a fraction of
% half the sampling rate.
Wn=2*[colo cohi]/Fs;

if Wn(2)>=1
  Wn(2)=0.99;
  warning('Frequencies adjusted to keep within the Nyquist rate')
end

if diff(Wn)<0.01
  warning(sprintf('%s\n%s\n%s\n%s',...
  ['Situations that seem to require an exceptionally narrow band'],...
  ['filter can be handled more reliably by decimation, filtering'],...
  ['with a filter of more  moderate band width, and interpolation'],...
  ['to the original sampling rate. (SAC manual)']))
end

[B,A]=feval(tipe,npol,Wn);

[H,F]=freqz(B,A,512,Fs);
HABS2=abs(H + eps).^2;

xf=filter(B,A,detrend(x(:)));

if npas==2
  xf=flipud(filter(B,A,detrend(flipud(xf(:)))));  
  HABS2=HABS2.^2;
end

warning off
EPB=bpmin(decibel(HABS2),F,3);
warning on

co=[colo cohi];
