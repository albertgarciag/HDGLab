function dv = stokes_2D_Ex2_analyticalVelocityD(X,problemParams,iMat) 
% Couette flow in a ring

% Points
x = X(:,1);
y = X(:,2);

% Parameters
a = -5/24;
b = 5/24;
r2 = x.^2+y.^2;
r4 = r2.^2;

% Gradient of the velocity
v1x = 2*b*x.*y./r4;
v1y = - a - b./r2 + 2*b*y.^2./r4;
v2x = a + b./r2 - 2*b*x.^2./r4;
v2y = -2*b*x.*y./r4;

dv = [v1x, v1y, v2x, v2y];