function dp = stokes_2D_Ex4_analyticalPressureD(X,problemParams,iMat) 

% Data
R = 1;
V = 1;

% Points
nOfPoints = size(X,1);

% Gradient
dp = zeros(nOfPoints,2);
dp(:,1) = -4*V/R^2;
