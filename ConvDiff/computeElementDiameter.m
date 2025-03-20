function [hMax,hElem] = computeElementDiameter(mesh)

nsd = size(mesh.X,2);
hMax = 0;
hElem = zeros(1,mesh.nOfElements);

if nsd == 2
    nOfElemVertices = 3;
elseif nsd == 3
    nOfElemVertices = 4;
end

for iElem = 1:mesh.nOfElements
    Te = mesh.indexT(iElem,1):mesh.indexT(iElem,2); 

    Tv = Te(1:nOfElemVertices);
    Xv = mesh.X(Tv,:);        

    for i = 1:nOfElemVertices-1
        Xi = Xv(i,:);
        for j = i+1:nOfElemVertices
            Xj = Xv(j,:);
            v = Xj - Xi;
            ds = sqrt(v*v');
            hMax = max(hMax, ds);
            hElem(iElem) = max(hElem(iElem),ds);
        end
    end
end