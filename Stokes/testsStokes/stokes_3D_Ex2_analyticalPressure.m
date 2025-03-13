function p = stokes_3D_Ex2_analyticalPressure(X,problemParams,iMat) 
% External flow past a sphere

% Data
R = 1;
Uinf = 1;
pinf = 0;
nu = problemParams.viscosity(iMat);

% Points
x = X(:,1);
y = X(:,2);
z = X(:,3);

r2 = x.^2 + y.^2 + z.^2;

% Pressure 
p = pinf - (3/2)*Uinf*nu*R*x./(r2.^(3/2));