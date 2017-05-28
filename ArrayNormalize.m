function y = ArrayNormalize(x, Dims, PrctLow, PrctHi)
% y = ArrayNormalize(x, Dims, PrctLow, PrctHi)
% 
% normalizes a multidimensional array according to percentiles.
%
% Dims says what dimensions to treat as independent. 
% for example if x was 3D and Dims was [] (the default), it would normalize 
% to the global max. If Dims was [3], it would normalize each (x,y) plane 
% individually. If Dims was [2,3], it would normalize each (x) line 
% individually. Etc.
%
% The data is linearly scaled so that the PrctLow'th percentile becomes 0 
% and the PrctHi'th percentile becomes 1. Defaults are 0 and 100.
% You might want to combine with clip()

if nargin<2; Dims = []; end
if nargin<3; PrctLow = 0; end
if nargin<4; PrctHi = 100; end
nDims = length(Dims);

Size = size(x);

SqueezeThese = setdiff(1:length(Size), Dims);

NewOrder = [Dims(:)' , SqueezeThese(:)'];
NewSize = [Size(Dims), prod(Size(SqueezeThese))];
if length(NewSize)==1
    NewSize = [NewSize, 1]; % I REALLY REALLY HATE MATLAB. NO, REALLY.
end

xp = permute(x, NewOrder);
xr = reshape(xp, NewSize);

Min = prctile(xr, PrctLow, nDims+1);
Max = prctile(xr, PrctHi, nDims+1);

rMin = ipermute(Min, NewOrder);
rMax = ipermute(Max, NewOrder);


y = bsxfun(@rdivide,bsxfun(@minus, double(x), double(rMin)), double(rMax-rMin));
