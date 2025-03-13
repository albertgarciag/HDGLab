function nodes = getNodesFace(nFace,nPerm,refElem)

nodes = refElem.face(nFace).nodesPerm(nPerm+1 ,:);