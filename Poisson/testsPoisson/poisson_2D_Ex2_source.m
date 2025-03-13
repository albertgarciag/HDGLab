function s = poisson_2D_Ex2_source(X, problemParams, iMat)

du2 = poisson_2D_Ex2_analyticalD2(X,problemParams,iMat);
k = problemParams.conductivity(iMat);
s = -k*( du2(:,1) + du2(:,4) );