function s = stokes_2D_Ex1_source(X, problemParams, iMat) 
% Wang flow in a square domain

dp  = stokes_2D_Ex1_analyticalPressureD(X,problemParams,iMat);
dv2 = stokes_2D_Ex1_analyticalVelocityD2(X,problemParams,iMat);
nu = problemParams.viscosity(iMat);

v1xx = dv2(:,1);
v1yy = dv2(:,4);
v2xx = dv2(:,5);
v2yy = dv2(:,8);

divGrad = [v1xx+v1yy, v2xx+v2yy];

s = dp-nu*divGrad;