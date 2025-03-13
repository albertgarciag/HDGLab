function dv = stokes_3D_Ex1_analyticalVelocityD(X,problemParams,iMat) 
% Synthetic solution in a unit cube

% Parameters
a = 1;
b = 0.5;

% Points
x = X(:,1);
y = X(:,2);
z = X(:,3);

% Gradient of the velocity
nOfPoints = numel(x);
dv = zeros(nOfPoints,9);

dv(:,1) = (z-y).*cos(x-b);
dv(:,2) = -sin(x-b);
dv(:,3) = sin(x-b);
dv(:,4) = y.*(z-0.5*y).*sin(x-b) - y.*cos(z-b);
dv(:,5) = -(z-y).*cos(x-b) - (x-y).*cos(z-b);
dv(:,6) = y.*(x-0.5*y).*sin(z-b) - y.*cos(x-b);
dv(:,7) = sin(z-b);
dv(:,8) = -sin(z-b);
dv(:,9) = (x-y).*cos(z-b);