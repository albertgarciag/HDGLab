function t = stokes_2D_Ex4_neumann(X, n, problemParams,iMat)

dv = stokes_2D_Ex4_analyticalVelocityD(X,problemParams,iMat);
p  = stokes_2D_Ex4_analyticalPressure(X,problemParams,iMat);
nu = problemParams.viscosity(iMat);

v1x = dv(:,1);
v1y = dv(:,2);
v2x = dv(:,3);
v2y = dv(:,4);

gradN1 = v1x.*n(:,1) + v1y.*n(:,2);
gradN2 = v2x.*n(:,1) + v2y.*n(:,2);
                                            
t = -n.*p+nu*[gradN1,gradN2];
