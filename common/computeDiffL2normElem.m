function [eL2, uL2] = computeDiffL2normElem(NuAtGaussRef,ue,refElemRef,XeRef,ueRef)

% Information of the reference element (reference solution)
NRef    = refElemRef.shapeFunctions(:,:,1)';

% Jacobian (reference solution)
if refElemRef.nsd==2
    [detJRef, ~] = isoparametricInvJ2D(XeRef, refElemRef.shapeFunctions(:,:,2), refElemRef.shapeFunctions(:,:,3));
elseif refElemRef.nsd==3
    [detJRef, ~] = isoparametricInvJ3D(XeRef, refElemRef.shapeFunctions(:,:,2), refElemRef.shapeFunctions(:,:,3), refElemRef.shapeFunctions(:,:,4));
end
weightsRef = refElemRef.gaussWeights.*detJRef;

% Information of the reference element (numerical solution)
N    = NuAtGaussRef(:,:,1)';

% Numerical and reference solutions at integration points
uG = N'*ue;
uGRef = NRef'*ueRef;

% Compute elemental contribution
eL2 = weightsRef'*(uG - uGRef).^2;
uL2 = weightsRef'*uGRef.^2;