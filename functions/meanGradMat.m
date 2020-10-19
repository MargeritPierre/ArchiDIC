function [D1,D2] = meanGradMat(mesh)
% Return sparse matrices D1 & D2 so that for any function f(X)
% so that df_dX1 = D1*[f] and df_dX2 = D2*[f]
% where f [nNodes 1] is the value of the function f at the mesh nodes
% and df_dX1 [nElems 1] is the mean gradient of the function f over each mesh element
%
% The mean gradient is computed using the Green theorem:
% <df_dx> = (1/A)*int(df_dx.da) = (1/A)*int(outer(f,n).dl)
% with n the outgoing normal of the contour
%
% The mesh has to be 2D and can contain either triangles OR quadrangles

[nElems,nMaxNodesByElem] = size(mesh.Elems) ;
[nNodes,nCoord] = size(mesh.Nodes) ;

% Element nodes coordinates
xe = reshape(mesh.Nodes(mesh.Elems(:),:),[nElems nMaxNodesByElem nCoord])  ;

% Element centroids
Ce = mean(xe,2,'omitnan') ; % [nElems 1 2]

% Element areas
areas = polyarea(xe(:,:,1),xe(:,:,2),2) ;

% Check element sorting (should be counter-clockwise)
xr = xe-Ce ; % relative coordinates w.r.t centroids [nElems nMaxNodesByElem 2]
xp = xr(:,:,1) + 1i*xr(:,:,2) ; % polar reprezentation [nElems nMaxNodesByElem 1]
da = angle(xp./circshift(xp,1,2)) ; % angle differences [nElems nMaxNodesByElem 1]
isSorted = all(da>=0,2) ; % all angle diff. should be positive [nElems 1 1]

% Element edges for integration
edges = cat(3,mesh.Elems,circshift(mesh.Elems,-1,2)) ; % [nElems nMaxNodesByElem 2]
xed = reshape(mesh.Nodes(edges(:),:),[nElems nMaxNodesByElem 2 nCoord]) ;
edgVec = diff(xed,1,3) ; % edge vectors [nElems nMaxNodesByElem 1 nCoord]

% Outgoing normals: edge vector turned by -pi/2 (outgoing if counter-clockwise sort)
normals = cat(3,edgVec(:,:,2),-edgVec(:,:,1)) ; % [nElems nMaxNodesByElem 2]
normals = normals.*(isSorted-0.5) ; % shift badly-oriented normals ;


% Build the sparse matrices
ee = repmat((1:nElems)',[1 nMaxNodesByElem 2]) ; % element indices
normals = normals./areas ; % integration divided by element area
n1 = repmat(normals(:,:,1),[1 1 2]) ; % first normal component
n2 = repmat(normals(:,:,2),[1 1 2]) ; % second normal component
D1 = sparse(ee(:),edges(:),n1(:),nElems,nNodes) ;
D2 = sparse(ee(:),edges(:),n2(:),nElems,nNodes) ;

end



function test
%% TEST
mesh = [] ;
mesh.Nodes = rand(100,2) ;
mesh.Elems = delaunay(mesh.Nodes) ;

[D1,D2] = meanGradMat(mesh) ;

fcn = mesh.Nodes(:,2) ;
[D1*fcn D2*fcn]

end




