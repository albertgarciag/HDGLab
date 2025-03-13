function p = stokes_2D_Ex1_analyticalPressure(X,problemParams,iMat) 
% Wang flow in a square domain

% Points
nOfPoints = size(X,1);

% Pressure 
p = zeros(nOfPoints, 1);