function t= poisson_3D_Ex1_Neumann(X, n, problemParams, iMat)

% NOTE that the normal is taken as an input so
% this is only valid for polygonal boundaries or NEFEM (exact normal)
% otherwise n should be defined analytically

du = poisson_3D_Ex1_analyticalD(X,problemParams,iMat);
k = problemParams.conductivity(iMat);
t = k*( du(:,1).*n(:,1) + du(:,2).*n(:,2) + du(:,3).*n(:,3) );