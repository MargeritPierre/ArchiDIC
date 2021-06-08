%% THIS SCRIPTS SHOWS HOW TO BUILD A TRIANGULATED MESH 
% USING DISTMESH APPLIED ON A REFERENCE IMAGE

%% MASK CREATION GUI
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
    roiType = 'poly' ; % 'rect' or 'poly'
    roiPosition = [72 77 734 587] ; [1 1 nJ-1 nI-1] ; % [xmin ymin width height]
    refPtPosition = [435 364] ; [nJ nI]/2 ;
    grayLevelThreshold = 0.075 ;
% Display parameters
    maskColor = [1 0 0] ;
    maskAlpha = 0.5 ;
    uiHeight = 0.04 ;
    lblWidth = 0.3 ;
    uiMargin = 0.005 ;
% Create the figure
    [fig,ax] = initFigure() ;
    fig.Name = 'MASK CREATION' ;
    ax.OuterPosition(4) = 1-3*uiHeight ;
% Reference image 
    im = image(refImg) ;
% Mask image
    maskIm = image(repmat(reshape(maskColor,[1 1 3]),[nI nJ])) ;
% Region Of Interest
    switch roiType
        case 'rect'
            roi = drawrectangle(ax,'Position',roiPosition,'FaceAlpha',0) ;
        case 'poly'
            roi = drawpolygon(ax,'Position',roiPosition(1:2) + [0 0 ; 1 0 ; 1 1 ; 0 1].*roiPosition(3:4) ,'FaceAlpha',0) ;
    end
% Reference point for the background gray level
    refPt = drawpoint(ax,'Position',refPtPosition) ;
% Slider for the threshold value
    tLabel = uicontrol(fig,'style','text'...
                            ,'units','normalized'...
                            ,'position',[0 1-uiHeight lblWidth uiHeight] + uiMargin*[1 1 -2 -2]) ;
    tSlider = uicontrol(fig,'style','slider'...
                            ,'units','normalized'...
                            ,'position',[lblWidth 1-uiHeight 1-lblWidth uiHeight] + uiMargin*[1 1 -2 -2]...
                            ,'min',0,'max',1,'value',grayLevelThreshold) ;
% Slider for the filter kernel radius
    RLabel = uicontrol(fig,'style','text'...
                            ,'units','normalized'...
                            ,'position',[0 1-2*uiHeight lblWidth uiHeight] + uiMargin*[1 1 -2 -2]) ;
    RSlider = uicontrol(fig,'style','slider'...
                            ,'units','normalized'...
                            ,'position',[lblWidth 1-2*uiHeight 1-lblWidth uiHeight] + uiMargin*[1 1 -2 -2]...
                            ,'min',1,'max',20,'value',10) ;
% Slider for the filter kernel order
    ordLabel = uicontrol(fig,'style','text'...
                            ,'units','normalized'...
                            ,'position',[0 1-3*uiHeight lblWidth uiHeight] + uiMargin*[1 1 -2 -2]) ;
    ordSlider = uicontrol(fig,'style','slider'...
                            ,'units','normalized'...
                            ,'position',[lblWidth 1-3*uiHeight 1-lblWidth uiHeight] + uiMargin*[1 1 -2 -2]...
                            ,'min',0,'max',1,'value',0.58) ;
% UI Listeners
    uiMask = @()img2mask(refImg,refPt.Position,roi.Position,tSlider.Value,RSlider.Value,ordSlider.Value) ;
    callbackFcn = @()cellfun(@(c)c(),{...
                        @()set(maskIm,'alphadata',uiMask()*maskAlpha) ...
                        @()set(tLabel,'string',['Gray level threshold: ' num2str(tSlider.Value)]) ...
                        @()set(RLabel,'string',['Filter Radius: ' num2str(round(RSlider.Value))]) ...
                        @()set(ordLabel,'string',['Filter order: ' num2str(round(ordSlider.Value)*100) '%']) ...
                        },'UniformOutput',false) ;
    addlistener(roi,'MovingROI',@(src,evt)callbackFcn()) ;
    addlistener(refPt,'MovingROI',@(src,evt)callbackFcn()) ;
    addlistener(tSlider,'Value','PostSet',@(src,evt)callbackFcn()) ;
    addlistener(RSlider,'Value','PostSet',@(src,evt)callbackFcn()) ;
    addlistener(ordSlider,'Value','PostSet',@(src,evt)callbackFcn()) ;
% Display
    callbackFcn() ;
    drawnow ;

%% LEVEL SET FUNCTION DEFINITION
% Parameters
    margin = 50 ; 
% Retrieve the mask
    MASK = uiMask() ;
% Bounding box
    bbox = round(roi.Position) ;
    switch roiType
        case 'rect'
            bbox = bbox(1:2)+[0;1]*bbox(3:4) ;
        case 'poly'
            bbox = [min(bbox,[],1) ; max(bbox,[],1)] ;
    end
    bbox = bbox + [-1;1]*margin ;
    jj = bbox(1,1):bbox(2,1) ;
    ii = bbox(1,2):bbox(2,2) ;
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
    distIm = imagesc(jj,ii,DIST,'tag',tag,'AlphaData',0.75) ;
    colorbar ;
    drawnow
    
%% MESH CONSTRUCTION
% Mesh Density Function
    edgeLength = 15 ;
    fh = @(p)huniform(p)*edgeLength ;
% Interpolated distance function
    fd = @(p)bilinearInterp(DIST,p-bbox(1,:)+1) ;
% Initial edge length
    h0 = edgeLength ;
% Fixed points (for eventual features)
    pFix = [] ;
% Compute the distmesh
    [p,t] = distmesh2d_modif(fd,fh,h0,bbox,pFix) ;
    
%% SAVE THE MESH
    Nodes = p ;
    Elems = t ;
    save(['data' filesep 'distMesh.mat'],'Nodes','Elems') ;
    

        
        