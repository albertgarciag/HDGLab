function dv = stokes_2D_Ex4_analyticalVelocityD(X, problemParams,iMat) 

% Data
R = 1;
V = 1;

% Points
nOfPoints = size(X,1);
r = X(:,2);

% Gradient of the velocity
v1z = zeros(nOfPoints,1);
v1r = -2*V*r/R^2;
v2z = zeros(nOfPoints,1);
v2r = zeros(nOfPoints,1);

dv = [v1z, v1r, v2z, v2r];