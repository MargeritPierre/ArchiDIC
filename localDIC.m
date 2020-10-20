%% LOCAL DIC FOLLOWING NODE FEATURES
    close all
    clearvars

%% IMPORT THE IMAGES
% Browse the "images" folder
    imgFiles = dir(['images' filesep '*.png']) ;
    nIMG = numel(imgFiles) ;
% Import images
    IMG = {} ;
    wtbr = waitbar(0,'Loading Images') ;
    for ii = 1:nIMG
        imgName = [imgFiles(ii).folder filesep imgFiles(ii).name] ;
        img = imread(imgName) ;
        if size(img,3)>1 ; img = sum(img/size(img,3),3,'native') ; end
        IMG{end+1} = img ;
        wtbr = waitbar(ii/nIMG,wtbr,['Loading Images (' num2str(ii) '/' num2str(nIMG) ')']) ;
    end
    delete(wtbr) ;
    drawnow ;
% Informations
    nIMG = numel(IMG) ;
    [nI,nJ] = size(IMG{1}) ;
% Image normalization function
    normalize = @(img)double(img)*(1/max(getrangefromclass(img(1)))) ;
    process = normalize ;
    
%% APPLY GAUSSIAN FILTERING IF REQUIRED (improves convergence)
    sigma = 0.6 ; filtSz = 11 ;
    filt = exp(-linspace(-1,1,filtSz).^2/sigma^2) ;
    filt = filt(:).*filt(:)' ;
    filt = filt/sum(filt(:)) ;
    process = @(img)conv2(normalize(img),filt,'same') ;
    
%% BUILD THE MACROSCOPIC QUAD MESH
% Parameters
    refImg = process(IMG{1}) ;
    gridSz = [2 3] ; % number of nodes in each dimension of the grid
    gridCornerPositions = [96 206 ; 774 200 ; 775 533 ; 97 539] ; [1 1 ; nJ 1 ; nJ nI ; 1 nI] ;
% Build the grid of points
    XI = arrayfun(@linspace,gridSz*0,gridSz*0+1,gridSz,'UniformOutput',false) ;
    [XI{:}] = ndgrid(XI{:}) ;
    XI = reshape(cat(numel(XI)+1,XI{:}),[],numel(XI)) ;
    pointsPosition = @(cornerPos)[(1-XI(:,1)).*(1-XI(:,2)) (1-XI(:,1)).*XI(:,2) XI(:,1).*XI(:,2) XI(:,1).*(1-XI(:,2))]*cornerPos ;
% Quad mesh
    p1 = (1:gridSz(1)-1)' + gridSz(1)*(0:gridSz(2)-2) ;
    elems = p1(:) + [0 gridSz(1) gridSz(1)+1 1] ;
% Display
    fig = clf('reset') ;
    ax = axes('Outerposition',[0 0 1 1]) ;
        hold on
        axis tight
        axis equal
        im = image(repmat(refImg,[1 1 3])) ;
    pts = patch('vertices',pointsPosition(gridCornerPositions)...
                ,'faces',elems ...
                ,'facealpha',0 ...
                ,'marker','.','markersize',25 ...
                ,'linewidth',2 ...
                ) ;
    cornerPoly = drawpolygon(ax,'Position',gridCornerPositions,'facealpha',0,'linewidth',0.01) ;
    addlistener(cornerPoly,'MovingROI',@(src,evt)set(pts,'vertices',pointsPosition(cornerPoly.Position))) ;
    drawnow ;
    
%% RECORD THE MESH AND DEFINE THE SUBDOMAIN SIZE
% Parameters
    subDomSz = [1 1]*61 ; % in pixels, [X Y]
% Retrieve mesh information
    mesh = [] ;
    mesh.Elems = elems ;
    mesh.Nodes = pts.Vertices ;
    nElems = size(mesh.Elems,1) ;
    nNodes = size(mesh.Nodes,1) ;
% Display Subdomains
    tag = 'subdomains' ;
    delete(findobj(gcf,'tag',tag)) ;
    boxes = mesh.Nodes + 0.5*subDomSz.*reshape([-1 -1 ; 1 -1 ; 1 1 ; -1 1]',[1 2 4]) ;
    boxes = permute(boxes,[1 3 2]) ;
    pa = patch('Vertices',reshape(boxes,[],2) ...
                ,'faces',(1:size(boxes,1))'+(0:3)*size(boxes,1) ...
                ,'EdgeColor','r' ...
                ,'LineWidth',2 ...
                ,'facealpha',0 ...
                ,'tag',tag ...
                ) ;
    drawnow ;
    
%% APPLY LOCAL DIC
% !! WARNING !!: THIS LOCAL DIC IMPLEMENTATION UPDATES ALL SUBDOMAIN POSITIONS AT
% ONCE: AS A CONSEQUENCE, THE IDENTIFIED DISPLACEMENTS ARE NO LONGER INDEPENDENT IF
% THE SUBDOMAINS OVERLAP !!
% SHAPE FUNCTIONS
    ii = round(mesh.Nodes(:,2) - 0.5*subDomSz(2) + (0:subDomSz(2)-1)) ; % image row indices
    jj = round(mesh.Nodes(:,1) - 0.5*subDomSz(1) + (0:subDomSz(1)-1)) ; % image column indices
    pp = repmat(ii,[1 size(jj,2)]) + nI*repelem(jj-1,1,size(ii,2)) ; % image pixel indices
    ROI = sort(pp(:)) ;
    N = sparse(pp(:)',repmat(1:nNodes,[1 size(pp,2)]),1,nI*nJ,nNodes) ;
    N = N(ROI,:) ;
% Reference configuration
    ii = mod(ROI-1,nI)+1 ;
    jj = (ROI-ii)/nI ;
    X = [jj(:) ii(:)] ;
% GRADIENT OF THE REFERENCE IMAGE dG_dX
    % Derivation kernels (finite differences)
        d_dX1 = [1;1;1]*[1 0 -1]/2 ;
        d_dX2 = [1;0;-1]*[1 1 1]/2 ;
    % Image gradient dG_dX
        dG_dX1 = conv2(refImg,d_dX1,'same') ;
        dG_dX2 = conv2(refImg,d_dX2,'same') ;
% PRE-COMPUTE THE HESSIAN
% with U = [u1(:) u2(:)] the nodal displacements to identify
% the vector of parameters is [U] = U(:) = [u1(:) ; u2(:)]
    dG_du = [dG_dX1(ROI).*N dG_dX2(ROI).*N] ; % [nPixels 2*nNodes] 
    H = dG_du'*dG_du ;
% PERFORM LOCAL DIC
    U = zeros(nNodes,2,nIMG) ;
% Parameters
    G = refImg(ROI) ; % Reference image on ROI
    U(:,:,1) = 0 ; % zero initial displacement
    addPreviousUpdate = true ; % helps convergence ;)
    criterion = @(dU)max(abs(dU(:))) ; % convergence criterion
    delta = 1e-3 ; % convergence threshold
    maxIt = 200 ; % maximum number of iterations
    plotFreq = 0 ; % display update frequency
    profile on
% Display
    fig = clf ;
    ax = axes('Outerposition',[0 0 1 1]) ;
        hold on
        im = image(repmat(refImg,[1 1 3])) ;
        pa = patch('Vertices',reshape(boxes,[],2) ...
                    ,'Faces',(1:size(boxes,1))'+(0:3)*size(boxes,1) ...
                    ,'FaceAlpha',0 ...
                    ,'EdgeColor','r' ...
                    ,'LineWidth',2 ...
                    ) ;
        axis tight
        axis equal
% For each image...
    ttt = tic ;
    for ii = 2:nIMG
    % Retrieve the current image
        img = process(IMG{ii}) ;
    % Intialize the displacement with the previous image
        U(:,:,ii) = U(:,:,ii-1) ;
        if addPreviousUpdate && ii>=3 ; U(:,:,ii) = U(:,:,ii) + U(:,:,ii-1) - U(:,:,ii-2) ; end
        it = 0 ; iterate = true ;
    % Iterate
        while iterate 
        % Current configuration at pixels
            xx = X + N*U(:,:,ii) ;
        % current image interpolation
            g = bilinearInterp(img,xx) ;
        % image residual
            res = g(:)-G(:) ;
        % Jacobian
            j = dG_du'*res ;
        % Update
            dU = H \ j ;
            dU = reshape(dU,[nNodes 2]) ;
            U(:,:,ii) = U(:,:,ii) - dU ;
        % Iterate
            it = it + 1 ;
            iterate = criterion(dU)>delta && it<maxIt ;
        % Display
            if it<2 ; im.CData = repmat(img,[1 1 3]) ; end
            if ~iterate || toc(ttt)>1/plotFreq
                pa.Vertices = reshape(boxes + permute(U(:,:,ii),[1 3 2]),[],2) ;
                fig.Name = [...
                             'Image ' num2str(ii) '/' num2str(nIMG) ' : ' ...
                             'iteration ' num2str(it) ...
                             ' , criterion = ', num2str(criterion(dU),'%g') ...
                           ] ;
                ax.XLim = [min(xx(:,1)) max(xx(:,1))] + subDomSz(1)*[-1 1] ;
                ax.YLim = [min(xx(:,2)) max(xx(:,2))] + subDomSz(2)*[-1 1] ;
                drawnow ;
                ttt = tic ;
            end
        end
    end
    profile off
    
    
%% PLOT THE DATA
    % Derivation matrix B such that EPSILON = B*[U]
        [D1,D2] = meanGradMat(mesh) ;
        O = D1*0 ;
        B = [D1 O ; O D2 ; D2 D1] ;
    % Data Fields (add fields if required)
        DataFields = [] ;
        % Displacement
        DataFields.u1 = U(:,1,:) ;
        DataFields.u2 = U(:,2,:) ;
        DataFields.U = sqrt(sum(U.^2,2)) ;
        % Linearized Strains
        EPS = B*reshape(U,[2*nNodes nIMG]) ; % [nElems*3 nIMG]
        EPS = reshape(EPS,[nElems 3 nIMG]) ; % [Nodes 3 nIMG]
        DataFields.E11 = EPS(:,1,:) ;
        DataFields.E22 = EPS(:,2,:) ;
        DataFields.E12 = 0.5*EPS(:,3,:) ;
    % Figure
        uiHeight = 0.03 ;
        popupWidth = 0.25 ;
        fig = clf('reset') ;
        fig.Name = 'Local DIC: Results' ;
        ax = axes('Outerposition',[0 0 1 1]) ;
            hold on
            im = image(repmat(refImg,[1 1 3])) ;
            pa = patch('Vertices',mesh.Nodes,'Faces',mesh.Elems,'FaceColor','interp','FaceAlpha',0.5) ;
            axis tight
            axis equal
        colorbar(ax) ;
        imgSlider = uicontrol(fig,'style','slider'...
                                ,'units','normalized'...
                                ,'position',[0 1-uiHeight 1-popupWidth uiHeight]...
                                ,'min',1,'max',nIMG,'value',1) ;
        dataPopup = uicontrol(fig,'style','popupmenu'...
                                ,'units','normalized'...
                                ,'position',[1-popupWidth 1-uiHeight popupWidth uiHeight]...
                                ,'String',fieldnames(DataFields)) ;
        plotFcn = @(ii,data)cellfun(@(fcn)fcn(ii,data),{...
                    @(ii,data)set(pa,'vertices',mesh.Nodes + U(:,:,ii),'FaceVertexCData',data) ...
                    @(ii,data)set(pa,'facecolor',char((size(data,1)==nNodes)*'interp' + (size(data,1)==nElems)*'flat  ')) ...
                    @(ii,data)set(im,'cdata',repmat(IMG{ii},[1 1 3]))' ...
                    @(ii,data)set(ax,'clim',[min(data(:))-eps max(data(:))+eps])' ...
                    },'UniformOutput',false) ;
        callbackFcn = @()plotFcn(round(imgSlider.Value),DataFields.(dataPopup.String{dataPopup.Value})(:,:,round(imgSlider.Value))) ;
        addlistener(imgSlider,'Value','PostSet',@(src,evt)callbackFcn()) ;
        dataPopup.Callback = @(src,evt)callbackFcn() ;
        callbackFcn() ;
        
%% SAVE THE RESULTS
    save(['data' filesep 'localDIC.mat'],'mesh','DataFields') ;







