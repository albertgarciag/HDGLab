function s = convdiff_2D_Ex3_source(X, problemParams, iMat)

a = convdiff_2D_Ex3_convection(X,problemParams,iMat);
du = convdiff_2D_Ex3_analyticalD(X,problemParams,iMat);
du2 = convdiff_2D_Ex3_analyticalD2(X,problemParams,iMat);
k = problemParams.conductivity(iMat);
s = a(:,1).*du(:,1) + a(:,2).*du(:,2)  -k*( du2(:,1) + du2(:,4) );