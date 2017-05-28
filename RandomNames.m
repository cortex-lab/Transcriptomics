function MyNames = RandomNames(n, NameFile);
% MyNames = RandomNames(n, NameFile);
% returns n randomly selected names sampled at random from a list of 2013
% baby names released by the US Social Security administration
%
% NameFile should point to a .mat file containing a single cell array Names

u = userpath;
if nargin<2
    if u(end)==';'
        NameFile = [u(1:end-1) '\BabyNames.mat'];
    else
        NameFile = [u '\BabyNames.mat'];
    end
end

load(NameFile)

p = randperm(length(Names));

MyNames = Names(p(1:n));