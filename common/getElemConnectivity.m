function Te = getElemConnectivity(mesh, iElem)

Te = mesh.indexT(iElem,1):mesh.indexT(iElem,2);