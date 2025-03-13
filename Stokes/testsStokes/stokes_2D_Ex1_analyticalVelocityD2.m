function dv2 = stokes_2D_Ex1_analyticalVelocityD2(X,problemParams,iMat) 
% Wang flow in a square domain

% Parameters
lambda = 10;
b = 1;

% Points
x = X(:,1);
y = X(:,2);

% Second derivatives of the velocity
v1xx = b*lambda^3*exp(-lambda*y).*cos(lambda*x);
v1xy = -b*lambda^3*exp(-lambda*y).*sin(lambda*x);
v1yx = v1xy;
v1yy = -v1xx;
v2xx = v1xy;
v2xy = v1yy;
v2yx = v2xy;
v2yy = -v1xy;

dv2 = [v1xx, v1xy, v1yx, v1yy, v2xx, v2xy, v2yx, v2yy];