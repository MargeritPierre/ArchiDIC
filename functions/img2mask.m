function [fullMask,img] = img2mask(img,refPt,ROI,t,R,ord)
%IMAGEMASK Create a binary mask from an image
% img: [nI nJ nC] image to process
% refPt: [x y] for background gray level substraction
% ROI: [xmin ymin width height] image cropping (or [x(:) y(:)] for polygonal ROI)
% t: scalar in [0->1] gray level threshold value
% R: scalar, statistical order filter radius
% ord: scalar in [0->1], statistical order threshold

% Convert to double
    if ~isa(img,'double')
        img = double(img)/max(getrangefromclass(img(1))) ;
    end
    
% Image infos
    [nI,nJ,~] = size(img) ;
    
% ROI info
    if size(ROI,1)==1 % rectangle [xmin ymin width height]
        ROI = ROI(1:2) + [0 0 ; 1 0 ; 1 1 ; 0 1].*ROI(3:4) ; % convert to polygon type
    else % polygon [x(:) y(:)]
    end

% Cropping indices
    ROI = round(ROI) ; % integer pixel coordinates
    bbox = [min(ROI,[],1) ; max(ROI,[],1)] ; % bounding box

% Background substraction
    refPt = round(refPt) ;
    refValue = img(refPt(2),refPt(1),:) ;
    img = abs(img-refValue) ;
    
% Cropping
    m = ceil(R/2) ; % margin for filter radius
    ii = max(bbox(1,2)-m,1):min(bbox(2,2)+m,nI) ;
    jj = max(bbox(1,1)-m,1):min(bbox(2,1)+m,nJ) ;
    img = img(ii,jj,:) ;
    
% Thresholding
    mask = mean(img,3)>t ;
    
% Statistical Order Filtering
    kernel = ( ((-R:R).^2)' + (-R:R).^2 )<=R.^2 ; % round disk binary kernel ;
    ord = round(ord.*sum(kernel(:))) ;
    mask = ordfilt2(mask,ord,kernel,'symmetric') ;
    
% ROI shape indices
    [JJ,II] = meshgrid(jj,ii) ;
    roi = inpolygon(JJ(:),II(:),ROI(:,1),ROI(:,2)) ;
    mask(~roi) = 0 ;
    
% Full mask
    fullMask = zeros([nI nJ],class(mask)) ;
    fullMask(ii,jj) = mask ;

end

