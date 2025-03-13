function s = stokes_3D_Ex1_source(X, problemParams, iMat) 
% Synthetic solution in a unit cube

dp  = stokes_3D_Ex1_analyticalPressureD(X,problemParams,iMat);
dv2 = stokes_3D_Ex1_analyticalVelocityD2(X,problemParams,iMat);
nu = problemParams.viscosity(iMat);

v1xx = dv2(:,1);
v1yy = dv2(:,5);
v1zz = dv2(:,9);
v2xx = dv2(:,10);
v2yy = dv2(:,14);
v2zz = dv2(:,18);
v3xx = dv2(:,19);
v3yy = dv2(:,23);
v3zz = dv2(:,27);

divGrad = [v1xx+v1yy+v1zz, v2xx+v2yy+v2zz, v3xx+v3yy+v3zz];

s = dp-nu*divGrad;