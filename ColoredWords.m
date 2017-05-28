function h = ColoredWords(Words, Colors, nColumns, ColorMap, Order)
% h = ColoredWords(Words, Colors, Columns, ColorMap, Order);
% 
% Makes a multicolor representation of some words in 2d
% Words: cell array of strings
% Colors: a number 0 to 1 for each word, saying what color it should be
% Columns: how many columns to plot (default ceil(sqrt(nWords)))
% ColorMap: default is upside-down hot
%
% Order is what order to plot them in: with the highest number on top. 
% Default is to sort so the largest value goes on top
%
% returns a vector of handles to the words

nWords = length(Words);

if nWords ==0
    h = [];
    return;
end

if nargin<3 | isempty(nColumns)
    nColumns = ceil(sqrt(nWords));
end

if nargin<4 | isempty(ColorMap);
    ColorMap = flipud(hot(200));
end

if nargin<5 | isempty(Colors)
    [~, Order] = sort(Colors);
end

if max(Colors)>1 || min(Colors)<0
    error('Colors should be numbers between 0 and 1');
end

Colors(~isfinite(Colors)) = .5;



nColors = size(ColorMap,1);

h = zeros(length(Order),1);

for i=Order(:)';
    h(i) = text(.5+mod(i-1,nColumns), -.5-floor((i-1)/nColumns), Words{i});
    set(h(i), 'Color', ColorMap(1+floor(Colors(i)*(nColors-1)),:));
    
end

axis([0, nColumns+1, -1-ceil(nWords/nColumns),0.5]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
set(gca, 'xticklabel', '')
set(gca, 'yticklabel', '')