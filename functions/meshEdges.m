function [edges,ele2edg] = meshEdges(mesh)
% Unique list of edges of the mesh
% edges [nEdges 2] node indices of each edge
% ele2edg [nEdges nElems] element->edge connectivity

[nElems,nMaxNodesByElem] = size(mesh.Elems) ;

% All mesh edges
edges = cat(3,mesh.Elems,circshift(mesh.Elems,-1,2)) ; % [nElems nMaxNodesByElem 2]
ee = repmat((1:nElems)',[1 nMaxNodesByElem 1]) ; % corresponding element indices

% Unique edges
edges = reshape(edges,[nElems*nMaxNodesByElem 2]) ; % [nElems*nMaxNodesByElem 2]
[edges,~,iedg] = unique(sort(edges,2),'rows') ;

% Edge connectivity
if nargout>1 
    ele2edg = sparse(iedg(:),ee(:),1,size(edges,1),nElems) ;
end

end

function test
%%

mesh.Nodes = rand(10,2) ;
mesh.Elems = delaunay(mesh.Nodes) ;

[edges,ele2edg] = meshEdges(mesh)


end

