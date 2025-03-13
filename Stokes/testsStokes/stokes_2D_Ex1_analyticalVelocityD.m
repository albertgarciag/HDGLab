function dv = stokes_2D_Ex1_analyticalVelocityD(X,problemParams,iMat) 
% Wang flow in a square domain

% Parameters
lambda = 10;
a = 1;
b = 1;

% Points
x = X(:,1);
y = X(:,2);

% Gradient of the velocity
v1x = b*lambda^2*exp(-lambda*y).*sin(lambda*x);
v1y = 2*a + b*lambda^2*exp(-lambda*y).*cos(lambda*x);
v2x = b*lambda^2*exp(-lambda*y).*cos(lambda*x);
v2y = -v1x;

dv = [v1x, v1y, v2x, v2y];