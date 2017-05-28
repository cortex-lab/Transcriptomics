function cm = RedWhiteBlue(n, gamma)
% RedWhiteBlue(n, gamma)
% colormap going from red to white to blue
% n is number of distinct colors
% gamma is a gamma correction term (use to emphasize small values)

if nargin<1; n = 100; end
if nargin<2; gamma = 0.6; end

cm = ([n*ones(1,n), n:-1:0 ; ...
      0:n, n-1:-1:0; ...
      0:n, ones(1,n)*n]' / n).^gamma;
% colormap(cm);
