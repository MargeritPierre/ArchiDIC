%% THIS SCRIPTS SHOWS HOW TO BUILD A TRIANGULATED MESH 
% USING DISTMESH APPLIED ON A REFERENCE IMAGE

%% MASK CREATION
% Clear the workspace
    close all
    clearvars
% Reference image
    imgFiles = dir(['images' filesep '*.png']) ;
    refImg = [imgFiles(1).folder filesep imgFiles(1).name] ;
    refImg = imread(refImg) ;
    refImg = double(refImg)/max(getrangefromclass(refImg)) ;
    [nI,nJ,nC] = size(refImg) ;
    if nC<3 ; refImg = repmat(refImg(:,:,1),[1 1 3]) ; end
% Parameters
    roiPosition = [72 77 734 587] ; [1 1 nJ-1 nI-1] ;
    refPtPosition = [435 364] ; [nJ nI]/2 ;
    grayLevelThreshold = 0.075 ;
    maskColor = [1 0 0] ;
    maskAlpha = 0.5 ;
    uiHeight = 0.04 ;
    uiMargin = 0.005 ;
% Create the figure
    fig = figure ;
        fig.NumberTitle = 'off' ;
        fig.Name = 'MASK CREATION: select ROI, substract background and apply threshold' ;
        fig.MenuBar = 'none' ; 
        fig.ToolBar = 'none' ;
    ax = axes('Outerposition',[0 0 1 1]) ;
        hold on ;
        axis(ax,'equal')
        axis(ax,'tight')
    im = image(refImg) ;
% Rectangular Region Of Interest
    roi = drawrectangle(ax,'Position',roiPosition,'FaceAlpha',0) ;
% Substract the background gray level from the image
    % Function
        substractBackground = @(pt)abs(refImg-refImg(round(pt(2)),round(pt(1)),:)) ;
    % Reference point for the background gray level
        refPt = drawpoint(ax,'Position',refPtPosition) ;
    % Listener (play with the point to see the image changing)
        set(im,'cdata',substractBackground(refPt.Position))
        addlistener(refPt,'MovingROI',@(src,evt)set(im,'cdata',substractBackground(src.Position))) ;
% MASK: Image thresholding
    % Function
        maskFun = @(t)all(im.CData>t,3) ;
    % slider
        tSlider = uicontrol(fig,'style','slider','units','normalized','position',[0 1-uiHeight 1 uiHeight] + uiMargin*[1 1 -2 -2],'min',0,'max',1,'value',grayLevelThreshold) ;
    % graphical representation of the mask
        maskIm = image(repmat(reshape(maskColor,[1 1 3]),[nI nJ])) ;
        maskIm.AlphaData = maskFun(tSlider.Value)*maskAlpha ;
    % listen to the slider
        addlistener(tSlider,'Value','PostSet',@(src,evt)set(maskIm,'alphadata',maskFun(tSlider.Value)*maskAlpha)) ;
    drawnow ;

%% MEDIAN FILTERING OF THE MASK
    Rf = 10 ; order = round(pi*Rf.^2/1.7) ;
% Round mask filter
    filtMask = ( ((-Rf:Rf).^2)' + (-Rf:Rf).^2 )<=Rf.^2 ;
    MASK = ordfilt2(maskFun(tSlider.Value),order,filtMask) ;
% Display
    maskIm.AlphaData = MASK*0.5 ;
    drawnow ;

%% LEVEL SET FUNCTION DEFINITION
% Bounding box
    margin = 50 ;
    bbox = round(roi.Position) ;
    bbox = bbox(1:2)+[0;1]*bbox(3:4) + [-1;1]*margin ;
    jj = bbox(1,1):bbox(2,1) ;
    ii = bbox(1,2):bbox(2,2) ;
    [JJ,II] = meshgrid(jj,ii) ;
% ROI with margins
    ROI = MASK(max(min(ii,nI),1),max(min(jj,nJ),1)) ;
    ROI(:,1:margin) = 0 ;
    ROI(:,end-margin+1:end) = 0 ;
    ROI(1:margin,:) = 0 ;
    ROI(end-margin+1:end,:) = 0 ;
% Level set image
    DIST = bwdist(ROI) - bwdist(~ROI) ;
    DIST = DIST + ROI.*0.5 - ~ROI*0.5 ;
% Display
    tag = 'lvlstimg' ;
    delete(findobj(fig,'tag',tag)) ;
    distIm = imagesc(jj,ii,DIST,'tag',tag) ;
    drawnow
    
%% MESH CONSTRUCTION
% Mesh Density Function
    edgeLength = 12 ;
    fh = @(p)huniform(p)*edgeLength ;
% Interpolated distance function
    fd = @(p)bilinearInterp(DIST,p-bbox(1,:)+1) ;
% Compute the distmesh
    profile on
    [p,t] = distmesh2d_modif(fd,fh,edgeLength,bbox,[]) ;
    profile off
    
%% SAVE THE MESH
    Nodes = p ;
    Elems = t ;
    save(['data' filesep 'distMesh.mat'],'Nodes','Elems') ;
    

        
        