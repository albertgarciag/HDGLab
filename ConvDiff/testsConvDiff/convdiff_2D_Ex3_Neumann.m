function t= convdiff_2D_Ex3_Neumann(X, n, problemParams, iMat)

% NOTE that the normal is taken as an input so
% this is only valid for polygonal boundaries or NEFEM (exact normal)
% otherwise t should be defined analytically

du = convdiff_2D_Ex3_analyticalD(X,problemParams,iMat);
k = problemParams.conductivity(iMat);
if problemParams.totalFluxNeumann == 1
    u = convdiff_2D_Ex3_analytical(X,problemParams,iMat);
    a = convdiff_2D_Ex3_convection(X,problemParams,iMat);
    du = -(a.*u - du);
end
t = k*( du(:,1).*n(:,1) + du(:,2).*n(:,2) );