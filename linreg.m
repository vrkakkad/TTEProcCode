%LINREG - estimate local slope
%
% dydi = linreg(y,klen,dim)
%
% LINREG returns the local slope of the matrix Y, using a kernel length of
% klen, in dimension dim. Dim defaults to the first non-singleton dimesion
% of Y. LINREG assumes even spacing between samples, so if you need dy/dx,
% you'll just need to divide the result by dx. 
%
% Version History:
% Created 3/17/14 Peter Hollender

function [dydi r2] = linreg(Y,klen,dim)
if ~exist('dim','var')
    dim = find(size(Y)>1,1,'first');
end
W = ~isnan(Y);
Y(~W) = 0;
k = flipud((-(klen-1)/2:(klen-1)/2)');
sz = 1:numel(size(Y));
sz(dim) = 1;
sz(1) = dim;
k = permute(k,sz);
k1 = ones(size(k));
Sw = convn(W,k1,'same');
Sx = convn(W,k,'same');
Sxy  = convn(Y,k,'same');
Sy = convn(Y,k1,'same');
Sxx = convn(W,k.^2,'same');
dydi = ((Sxy-(Sx.*Sy)./Sw)./(Sxx-(Sx.^2./Sw)));
if nargout>1
Syy = convn(Y.^2,k1,'same');
r2 = (((Sxy-(Sx.*Sy)./Sw)./sqrt((Sxx-(Sx.^2)./Sw).*(Syy./Sw-(Sy./Sw).^2))).^2)./Sw;
end