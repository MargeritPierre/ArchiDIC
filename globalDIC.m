%% GLOBAL DIC USING A P1 TRIANGULATION MESH
    close all
    clearvars

%% IMPORT THE MESH
    mesh = load(['data' filesep 'Mesh.mat']) ;
    nNodes = size(mesh.Nodes,1) ;
    nElems = size(mesh.Elems,1) ;

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

%% COMPUTE THE SHAPE FUNCTIONS OF THE MESH
% Evaluated at pixel coordinates
    refImg = normalize(IMG{1}) ;
    N = P1ShapeFunctions(mesh,refImg) ; % shape functions
    ROI = find(any(N,2)) ; % Region of Interest
    N = N(ROI,:) ; % cull out-of ROI pixels

%% EVALUATE THE GRADIENT OF THE REFERENCE IMAGE
% Derivation kernels (finite differences)
    d_dX1 = [1;1;1]*[1 0 -1]/2 ;
    d_dX2 = [1;0;-1]*[1 1 1]/2 ;
% Image gradient dG_dX
    dG_dX1 = conv2(refImg,d_dX1,'same') ;
    dG_dX2 = conv2(refImg,d_dX2,'same') ;

%% PRE-COMPUTE THE HESSIAN
% with U = [u1(:) u2(:)] the nodal displacements to identify
% the vector of parameters is [U] = U(:) = [u1(:) ; u2(:)]
    dG_du = [dG_dX1(ROI).*N dG_dX2(ROI).*N] ; % [nPixels 2*nNodes] 
    H = dG_du'*dG_du ;
    
%% STRAIN REGULARIZATION
% Gradient matrices D1 & D2
    [D1,D2] = meanGradMat(mesh) ;
% Strain matrix B so that [E11(:) ; E22(:) ; 2*E12(:)] = B*[U]
    O = sparse(nElems,nNodes) ; % matrix full of zeros
    B = [D1 O ; O D2 ; D2 D1] ;
% Edge differentiation matrix G so that df = G*[f] where f is defined on
% elements and df measures the difference between 2 element values 
    % element->edge connectivity
    [~,ele2edg] = meshEdges(mesh) ;
    % Keep interior edges only 
    boundaryEdges = sum(ele2edg,2)==1 ;
    ele2edg(boundaryEdges,:) = [] ;
    % Difference matrix
    [elem,edg] = find(ele2edg') ;
    val = repmat([1;-1],[numel(edg)/2 1]) ;
    G = sparse(edg,elem,val,size(ele2edg,1),nElems) ;
% Distance between element centroids
    C = reshape(mean(reshape(mesh.Nodes(mesh.Elems(:),:),[nElems 3 2]),2),[nElems 2]) ;
    dC = sqrt(sum((G*C).^2,2)) ;
% Strain Gap
    Gc = (1./dC(:)).*G ; % difference normalized by the centroid distance
    dB = blkdiag(Gc,Gc,Gc)*B ; % dB*[U] measures the normalized strain gap
% Associated Hessian
    Hr = dB'*dB ;
    
%% PERFORM GLOBAL DIC (WITHOUT REGULARIZATION)
    U = zeros(nNodes,2,nIMG) ;
% Parameters
    G = refImg(ROI) ; % Reference image on ROI
    U(:,:,1) = 0 ; % zero initial displacement
    addPreviousUpdate = true ; % helps convergence ;)
    beta = 1e5 ; % regularization coefficient
    criterion = @(dU)max(abs(dU(:))) ; % convergence criterion
    delta = 1e-3 ; % convergence threshold
    maxIt = 200 ; % maximum number of iterations
    plotFreq = Inf ; % display update frequency
    profile on
% Display
    fig = clf ;
    ax = axes('Outerposition',[0 0 1 1]) ;
        im = image(repmat(refImg,[1 1 3]),'interpolation','bilinear') ;
        pa = patch('Vertices',mesh.Nodes,'Faces',mesh.Elems,'FaceAlpha',0,'EdgeColor','r') ;
        axis tight
        axis equal
% For each image...
    ttt = tic ;
    for ii = 2:nIMG
    % Retrieve the current image
        img = normalize(IMG{ii}) ;
    % Intialize the displacement with the previous image
        U(:,:,ii) = U(:,:,ii-1) ;
        if addPreviousUpdate && ii>=3 ; U(:,:,ii) = U(:,:,ii) + U(:,:,ii-1) - U(:,:,ii-2) ; end
        it = 0 ; iterate = true ;
    % Iterate
        while iterate 
        % Current configuration
            x = mesh.Nodes + U(:,:,ii) ; % nodal values
            xx = N*x ; % projection on pixels
        % current image interpolation
            g = bilinearInterp(img,xx) ;
        % image residual
            res = g(:)-G(:) ;
        % Jacobian
            j = dG_du'*res ;
            jr = Hr*reshape((U(:,:,ii)-U(:,:,ii-1)),[],1) ;
        % Update
            dU = (H + beta*Hr) \ (j + beta*jr) ;
            dU = reshape(dU,[nNodes 2]) ;
            U(:,:,ii) = U(:,:,ii) - dU ;
        % Iterate
            it = it + 1 ;
            iterate = criterion(dU)>delta && it<maxIt ;
        % Display
            if it<2 ; im.CData = repmat(img,[1 1 3]) ; end
            if ~iterate || toc(ttt)>1/plotFreq
                pa.Vertices = mesh.Nodes + U(:,:,ii) ;
                fig.Name = [...
                             'Image ' num2str(ii) '/' num2str(nIMG) ' : ' ...
                             'iteration ' num2str(it) ...
                             ' , criterion = ', num2str(criterion(dU),'%g') ...
                           ] ;
                ax.XLim = [min(xx(:,1)) max(xx(:,1))] + 0.05*[-1 1]*range(xx(:,1)) ;
                ax.YLim = [min(xx(:,2)) max(xx(:,2))] + 0.05*[-1 1]*range(xx(:,2)) ;
                drawnow ;
                ttt = tic ;
            end
        end
    end
    profile off
    
    
%% PLOT THE DATA
    % Node relocalization matrix: mean over attached elements
        ele2nod = sparse(mesh.Elems(:)',repmat(1:nElems,[1 3]),1,nNodes,nElems) ;
        ele2nod = (1./sum(ele2nod,2)).*ele2nod ;
    % Data Fields (add fields if required)
        DataFields = [] ;
        % Displacement
        DataFields.u1 = U(:,1,:) ;
        DataFields.u2 = U(:,2,:) ;
        DataFields.U = sqrt(sum(U.^2,2)) ;
        % Linearized Strains
        EPS = B*reshape(U,[2*nNodes nIMG]) ; % [nElems*3 nIMG]
        EPS = ele2nod*reshape(EPS,[nElems 3*nIMG]) ; % [Nodes 3*nIMG]
        EPS = reshape(EPS,[nNodes 3 nIMG]) ; % [nNodes 3 nIMG]
        DataFields.E11 = EPS(:,1,:) ;
        DataFields.E22 = EPS(:,2,:) ;
        DataFields.E12 = 0.5*EPS(:,3,:) ;
    % Figure
        uiHeight = 0.03 ;
        popupWidth = 0.25 ;
        fig = clf('reset') ;
        fig.Name = 'Global DIC: Results' ;
        ax = axes('Outerposition',[0 0 1 1]) ;
            im = image(repmat(refImg,[1 1 3]),'interpolation','bilinear') ;
            pa = patch('Vertices',mesh.Nodes,'Faces',mesh.Elems,'FaceColor','interp') ;
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
                    @(ii,data)set(im,'cdata',repmat(IMG{ii},[1 1 3]))' ...
                    @(ii,data)set(ax,'clim',[min(data(:))-eps max(data(:))+eps])' ...
                    },'UniformOutput',false) ;
        callbackFcn = @()plotFcn(round(imgSlider.Value),DataFields.(dataPopup.String{dataPopup.Value})(:,:,round(imgSlider.Value))) ;
        addlistener(imgSlider,'Value','PostSet',@(src,evt)callbackFcn()) ;
        dataPopup.Callback = @(src,evt)callbackFcn() ;
        callbackFcn() ;







