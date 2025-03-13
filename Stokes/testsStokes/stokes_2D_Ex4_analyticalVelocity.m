function v = stokes_2D_Ex4_analyticalVelocity(X,problemParams,iMat) 

% Data
R = 1;
V = 1;

% Points
nOfPoints = size(X,1);
r = X(:,2);

% Velocity
v1 =  V*(1-r.^2/R^2);
v2 = zeros(nOfPoints,1);

v = [v1, v2];