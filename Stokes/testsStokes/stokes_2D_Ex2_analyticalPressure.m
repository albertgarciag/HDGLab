function p = stokes_2D_Ex2_analyticalPressure(X,problemParams,iMat) 
% Couette flow in a ring

% Points
nOfPoints = size(X,1);

% Pressure 
p = ones(nOfPoints, 1);