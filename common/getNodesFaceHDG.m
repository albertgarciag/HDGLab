function nodes = getNodesFaceHDG(nFace,nPerm,refElem)

nodes = refElem.face(nFace).nodesPermHDG(nPerm+1 ,:);