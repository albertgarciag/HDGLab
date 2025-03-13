function dv2 = stokes_3D_Ex1_analyticalVelocityD2(X,problemParams,iMat) 
% Synthetic solution in a unit cube

% Parameters
b = 0.5;

% Points
x = X(:,1);
y = X(:,2);
z = X(:,3);

% Second derivatives of the velocity
nOfPoints = numel(x);
dv2 = zeros(nOfPoints,27);

dv2(:,1) = -(z-y).*sin(x-b);
dv2(:,2) = -cos(x-b);
dv2(:,3) = cos(x-b);
dv2(:,4) = -cos(x-b);
dv2(:,5) = 0;
dv2(:,6) = 0;
dv2(:,7) = cos(x-b);
dv2(:,8) = 0;
dv2(:,9) = 0;
dv2(:,10) = y.*(z-0.5*y).*cos(x-b);
dv2(:,11) = (z-y).*sin(x-b) - cos(z-b);
dv2(:,12) = y.*(sin(x-b) + sin(z-b));
dv2(:,13) = (z-y).*sin(x-b) - cos(z-b);
dv2(:,14) = cos(x-b) + cos(z-b);
dv2(:,15) = (x-y).*sin(z-b) - cos(x-b);
dv2(:,16) = y.*(sin(x-b) + sin(z-b));
dv2(:,17) = (x-y).*sin(z-b) - cos(x-b);
dv2(:,18) = y.*(x-0.5*y).*cos(z-b);
dv2(:,19) = 0;
dv2(:,20) = 0;
dv2(:,21) = cos(z-b);
dv2(:,22) = 0;
dv2(:,23) = 0;
dv2(:,24) = -cos(z-b);
dv2(:,25) = cos(z-b);
dv2(:,26) = -cos(z-b);
dv2(:,27) = -(x-y).*sin(z-b);