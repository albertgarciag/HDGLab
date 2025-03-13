function mesh2 = buildIsoparametricHOmesh(mesh, refElem)


mesh2 = mesh;
mesh2.pElem = mesh.pElem+1;

% Initialisation
nOfX = 0;
mesh2.indexT = zeros(mesh2.nOfElements,2);
indexIni = 1;

for iElem = 1:mesh2.nOfElements
    pElem = mesh2.pElem(iElem);
    indexEnd = refElem(pElem).nOfNodes;
    nOfX = nOfX + (indexEnd-indexIni+1);
end

mesh2.X = zeros(nOfX, mesh2.nsd);

% Computation
mesh2.indexT = zeros(mesh2.nOfElements,2);
indexIni = 1;
for iElem = 1:mesh2.nOfElements
    pElem = mesh2.pElem(iElem);
    indexEnd = indexIni + refElem(pElem).nOfNodes - 1;
    
    mesh2.indexT(iElem,1) = indexIni;
    mesh2.indexT(iElem,2) = indexEnd;
    
    indexIniOld = mesh.indexT(iElem,1);
    indexEndOld = mesh.indexT(iElem,2);
    XeOld = mesh.X(indexIniOld:indexEndOld,:);
    
    mesh2.X(indexIni:indexEnd,:) = refElem(pElem-1).shapeFunctionsNodesPPp1*XeOld;
    
    indexIni = indexEnd + 1;
end

mesh2.nOfNodes = size(mesh2.X,1);

