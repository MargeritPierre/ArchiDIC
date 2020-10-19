function N = P1ShapeFunctions(mesh,refImg)
% Compute the value of P1 triangular element shape functions at pixel
% coordinates of the reference image
% reImg [nI nJ 1] gray level image
% mesh: structure with at least the 2 fields:
%   - Nodes [nNodes 2]: positions of the mesh nodes in the image
%   - Elems [nElems 3]: triangulation (node indices)
% output N: [nPixels nNodes] (HUGE but SPARSE matrix !)

%tic
%clc,clearvars -except mesh refImg

% Infos
nNodes = size(mesh.Nodes,1) ;
nElems = size(mesh.Elems,1) ;
[nI,nJ,nC] = size(refImg) ;

% 1) In which element bounding box is which pixel ?
% Element bounding boxes
    xe = mesh.Nodes(mesh.Elems(:),:) ; % [nElems*3 2]
    xe = reshape(xe,[nElems 3 2]) ; % [nElems 3 2] 
    xmin = min(floor(xe),[],2) ; % [nElems 1 2]
    xmax = max(ceil(xe),[],2) ; % [nElems 1 2]
% concatenate all pixels
    nPixInBBoxDim = xmax-xmin+1 ; % [nElems 1 2]
    nPixInBBox = prod(nPixInBBoxDim,3) ; % [nElems 1 1]
    nPixToTest = sum(nPixInBBox,1) ; % [1 1 1]
    ee = repelem(1:nElems,nPixInBBox) ; % element index
    nPixInPreviousBBoxes = [0 cumsum(nPixInBBox(1:end-1)')] ;
    pp = (0:nPixToTest-1) - repelem(nPixInPreviousBBoxes,nPixInBBox) ; % pixel index
    ii = mod(pp,repelem(nPixInBBoxDim(:,2)',nPixInBBox)) ; % bbox row index
    jj = (pp-ii)./repelem(nPixInBBoxDim(:,2)',nPixInBBox) ; % bbox column index
    ii = ii + xmin(ee,2)' ; % image row index
    jj = jj + xmin(ee,1)' ; % image column index
    
% 2) Compute the shape functions
% Relative pixel coordinates
    PP = [jj(:) ii(:)] ; % Pixel coordinates [nPixToTest 2]
    xe = permute(xe,[1 3 2]) ; % triangle nodes [nElems 2 3]
    PP = PP-xe(ee,:,1) ; % Relative pixel coordinates [nPixToTest 2]
% Local coordinates xi such that A.xi = PP so xi = inv(A).PP
    A = xe(:,:,[2 3])-xe(:,:,1) ; % edge vectors [nElems 2 2]
    detA = A(:,1,1).*A(:,2,2)-A(:,2,1).*A(:,1,2) ; % [nElems 1 1]
    invA = [A(:,2,2) -A(:,2,1) -A(:,1,2) A(:,1,1)]./detA ; % [nElems 4]
    invA = reshape(invA,[nElems 2 2]) ; % [nElems 2 2] ;
    xi = sum(invA(ee,:,:).*reshape(PP,[nPixToTest 1 2]),3) ; % [nPixToTest 2]
% Shape Functions N = [1-Xi1-Xi2 Xi1 Xi2]
    NN = [1-xi(:,1)-xi(:,2) xi(:,1) xi(:,2)] ;
    
% 3) keep only valid tests (so that shape functions are in [0 1])
    valid = all(NN>=0 & NN<=1,2) ;
    ii = ii(valid) ;
    jj = jj(valid) ;
    NN = NN(valid,:) ;
    ee = ee(valid) ;
    
% 3) Build the sparse matrix
% Pixel indices
    nPix = nI*nJ ;
    pp = ii+(jj-1)*nI ;
% Node indices
    nn = mesh.Elems(ee,:) ;
% Sparse
    N = sparse(repmat(pp(:),[3 1]),nn(:),NN(:),nPix,nNodes) ;
    
%toc
end
