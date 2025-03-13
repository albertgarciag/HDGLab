function v = stokes_2D_Ex1_analyticalVelocity(X,problemParams,iMat) 
% Wang flow in a square domain

% Parameters
lambda = 10;
a = 1;
b = 1;

% Points
x = X(:,1);
y = X(:,2);

% Velocity
v1 = 2*a*y-b*lambda*exp(-lambda*y).*cos(lambda*x);
v2 = b*exp(-lambda*y).*sin(lambda*x)*lambda;

v = [v1, v2];