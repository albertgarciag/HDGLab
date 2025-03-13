function dp = stokes_2D_Ex1_analyticalPressureD(X,problemParams,iMat) 
% Wang flow in a square domain

% Points
nOfPoints = size(X,1);

% Gradient
dp = zeros(nOfPoints,2);