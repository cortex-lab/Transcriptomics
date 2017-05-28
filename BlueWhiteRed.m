function cm = BlueWhiteRed(varargin)
% BlueWhiteRed(n, gamma)
% colormap going from blue to white to red
% n is number of distinct colors
% gamma is a gamma correction term (use to emphasize small values)

cm = flipud(RedWhiteBlue(varargin{:}));