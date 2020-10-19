function G = edgeDiffMat(mesh)
% Edge differentiation matrix G so that df = G*[f] 
% where f is defined on elements 
% and df measures the gap in f between two adjacent elements

[nElems,nMaxNodesByElem] = size(mesh.Elems) ;
[nNodes,nCoord] = size(mesh.Nodes) ;

% All mesh edges
edges = cat(3,mesh.Elems,circshift(mesh.Elems,-1,2)) ; % [nElems nMaxNodesByElem 2]
ee = repmat((1:nElems)',[1 nMaxNodesByElem 2]) ; % corresponding element indices

% Unique edges
edges = reshape(edges,[nElems*nMaxNodesByElem 2],2) ; % [nElems*nMaxNodesByElem 2]
[edges,iedg,~] = unique(edges,'rows') ; 
nEdges = size(edges,1) ;

% Cull boundary edges (edges with only one attached element
G = sparse(iedg(:),ee(:),1,nEdges,nElems) ;
nAttachedElems = sum(G,2) ;
G(nAttachedElems<2,:) = 0 ;






end

