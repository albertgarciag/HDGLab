function t= poisson_2D_Ex2_Neumann(X, n, problemParams, iMat)

% NOTE that the normal is taken as an input so
% this is only valid for polygonal boundaries or NEFEM (exact normal)
% otherwise t should be defined analytically

du = poisson_2D_Ex2_analyticalD(X,problemParams,iMat);
k = problemParams.conductivity(iMat);
t = k*( du(:,1).*n(:,1) + du(:,2).*n(:,2) );