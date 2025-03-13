function v = stokes_3D_Ex2_analyticalVelocity(X,problemParams,iMat) 
% External flow past a sphere

v = zeros(size(X));

% Analytical solution of the flow past a sphere with velocity along x
% direction

% Data
R = 1;
Uinf = 1;

% Points
x = X(:,1);
y = X(:,2);
z = X(:,3);

r2 = x.^2 + y.^2 + z.^2;

% Velocity along x direction
v1 = Uinf + (1/4)*Uinf*R^3.*(3*x.^2./r2-1)./(r2.^(3/2)) - (3/4)*Uinf*R*(x.^2./r2+1)./(r2.^(1/2));
v2 = (3/4)*Uinf*R*y.*x.*(R^2./r2-1)./(r2.^(3/2));
v3 = (3/4)*Uinf*R*z.*x.*(R^2./r2-1)./(r2.^(3/2));

v = [v1, v2, v3];