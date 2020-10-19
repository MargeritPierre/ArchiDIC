function val = bilinearInterp(img,xx,extrapVal)
% Bilinear interpolation of an image at real-valued coordinates xx
% faster than gg = interp2(img,xx(:,1),xx(:,2),'linear') ;
% xx = [nValues 2] = [jj(:) ii(:)] ;
% with ii and jj resp. the (real-valued) row and column indices
% extrapVal: extrapolation value (default: NaN)

if nargin<3 ; extrapVal = NaN ; end

% Image informations
[nI,nJ,~] = size(img) ;

% Valid coordinates
valid = xx(:,1)<=nJ ...
        & xx(:,1)>=1 ...
        & xx(:,2)<=nI ...
        & xx(:,2)>=1 ;
    
% Dummy values
xx(~valid,:) = 1 ;

% Integer part of the cordinates
ji = floor(xx) ;
ji(ji(:,1)==nJ,1) = nJ-1 ;
ji(ji(:,2)==nI,2) = nI-1 ;

% Neightboring pixels
p1 = ji(:,2)+nI*(ji(:,1)-1) ;
p2 = p1 + nI ; 
p3 = p2+1 ; 
p4 = p1+1 ;

% residual coordinates
dx = xx-ji ;

% bilinear interpolation
val = img(p1).*(1-dx(:,1)).*(1-dx(:,2)) ...
    + img(p2).*dx(:,1).*(1-dx(:,2)) ...
    + img(p3).*dx(:,1).*dx(:,2) ...
    + img(p4).*(1-dx(:,1)).*dx(:,2) ;

% Non-valid values
val(~valid) = extrapVal ;

end

