function [fullMask,img] = img2mask(img,refPt,ROI,t,R,ord)
%IMAGEMASK Create a binary mask from an image
% img: [nI nJ nC] image to process
% refPt: [x y] for background gray level substraction
% ROI: [xmin ymin width height] image cropping
% t: scalar in [0->1] gray level threshold value
% R: scalar, statistical order filter radius
% ord: scalar in [0->1], statistical order threshold

% Convert to double
    if ~isa(img,'double')
        img = double(img)/max(getrangefromclass(img(1))) ;
    end
    
% Image infos
    [nI,nJ,nC] = size(img) ;

% Cropping indices
    ROI = round(ROI(:)') ;
    ROI = repmat(ROI(1:2),[1 2]) + kron([0 1],ROI(3:4)) ; %[xmin ymin ; xmax ymax]
    ii = ROI(2):ROI(4) ;
    jj = ROI(1):ROI(3) ;

% Background substraction
    refPt = round(refPt) ;
    refValue = img(refPt(2),refPt(1),:) ;
    img = abs(img-refValue) ;
    
% Thresholding
    mask = mean(img(ii,jj,:),3)>t ;
    
% Statistical Order Filtering
    R = round(R) ;
    kernel = ( ((-R:R).^2)' + (-R:R).^2 )<=R.^2 ; % round disk binary kernel ;
    ord = round(ord.*sum(kernel(:))) ;
    mask = ordfilt2(mask,ord,kernel,'symmetric') ;
    
% Full mask
    fullMask = zeros(size(img(:,:,1)),class(mask)) ;
    fullMask(ii,jj) = mask ;

end

