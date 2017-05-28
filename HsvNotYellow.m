function cm = HsvNotYellow(varargin)
% HsvNotYellow(n) makes a HSV colormap (see hsv) but darkens yellow colors
% so you can see them on a white background.

Thresh = .7;
Mult = 1;

% don't do 6 since the two ends are too similar
% if mod(varargin{1},5)==1
%     varargin{1}=varargin{1}+floor(varargin{1}/5);
% end

cm0 = hsv(varargin{:});
Yellowness = cm0(:,1)+cm0(:,2);
Divisor = 1 + Mult*max(Yellowness - Thresh, 0);
cm = bsxfun(@rdivide, cm0, Divisor);
