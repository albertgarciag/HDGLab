function t = stokes_3D_Ex1_neumann(X, n, problemParams,iMat)
% Synthetic solution in a unit cube

dv = stokes_3D_Ex1_analyticalVelocityD(X,problemParams,iMat);
p  = stokes_3D_Ex1_analyticalPressure(X,problemParams,iMat);
nu = problemParams.viscosity(iMat);

u1x = dv(:,1);
u1y = dv(:,2);
u1z = dv(:,3);
u2x = dv(:,4);
u2y = dv(:,5);
u2z = dv(:,6);
u3x = dv(:,7);
u3y = dv(:,8);
u3z = dv(:,9);

gradN1 = u1x.*n(:,1) + u1y.*n(:,2) + u1z.*n(:,3);
gradN2 = u2x.*n(:,1) + u2y.*n(:,2) + u2z.*n(:,3);
gradN3 = u3x.*n(:,1) + u3y.*n(:,2) + u3z.*n(:,3);

t = -n.*p+nu*[gradN1,gradN2,gradN3];