function dbscaled=decibel(data,ref)
% dbscaled=DECIBEL(data,ref)
%
% Scales data on a decibel scale,
% i.e. 10*log10(data/ref)
%
% INPUT:
%
% data       Can be a vector or a matrix
% ref        The reference value [default: max(data(:))]
%
% OUTPUT:
%
% dbscaled   Will retain the correct dimensions
% 
% Last modified by fjsimons-at-alum.mit.edu, 07/07/2008

[m,n]=size(data);

defval('ref',max(data(:)));
dbscaled=reshape(10*log10(data/ref),m,n);

% disp(sprintf('Using 10 log 10 (data / ref) with ref = %8.5f',ref))

% If you want to define power as pressure^2, say, then, dealing in
% pressure you'd need the prefactor to be 20. The current scale is good
% for units of power. dB SPL is 20log(p/0.0002 dynes/cm2) for comparisons
% with Sound Pressure Level.
