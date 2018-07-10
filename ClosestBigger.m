function Dist = ClosestBigger(d, v)
% [Dist] = ClosestBigger(d, v)
%
% given a square matrix d of distances, and a vector v of values,
% searches for the distance from each point to its nearest neighbor
% with a larger value.
%
% if v is a density, this does cluster analysis (see Rodriguez&Laio, Science 2014)
%
% Dist is distance to closest bigger (returns largest dist for top score)

n = size(d,1);

MaxDist = max(d(:));

d(1:(n+1):n^2) = MaxDist;

Dist = zeros(n,1);
Neighbor = zeros(n,1);

for i=1:n
    Bigger = find(v>v(i));
    if isempty(Bigger)
        Dist(i) = MaxDist;
    else
        Dist(i) = min(d(Bigger, i));
    end
end