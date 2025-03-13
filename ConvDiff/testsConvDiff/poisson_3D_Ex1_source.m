function s = poisson_3D_Ex1_source(X, problemParams, iMat)

du2 = poisson_3D_Ex1_analyticalD2(X,problemParams,iMat);
k = problemParams.conductivity(iMat);
s = -k*( du2(:,1) + du2(:,5) + du2(:,9) );