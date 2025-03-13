function v = stokes_3D_Ex1_analyticalVelocity(X,problemParams,iMat) 
% Synthetic solution in a unit cube

% Parameters
a = 1;
b = 0.5;

% Points
x = X(:,1);
y = X(:,2);
z = X(:,3);

% Velocity
v1 = b + (z-y).*sin(x-b);
v2 = a - y.*(z-0.5*y).*cos(x-b) - y.*(x-0.5*y).*cos(z-b);
v3 = b + (x-y).*sin(z-b);

v = [v1, v2, v3];