function v = stokes_2D_Ex2_analyticalVelocity(X,problemParams,iMat) 
% Couette flow in a ring

% Points
x = X(:,1);
y = X(:,2);

% Parameters
a = -5/24;
b = 5/24;
r2 = x.^2+y.^2;

% Velocity
v1 = -a*y - b*y./r2;
v2 = a*x + b*x./r2;

v = [v1, v2];