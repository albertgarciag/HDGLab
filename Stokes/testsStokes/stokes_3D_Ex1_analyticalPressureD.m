function dp = stokes_3D_Ex1_analyticalPressureD(X,problemParams,iMat) 
% Synthetic solution in a unit cube

% Points
x = X(:,1);
y = X(:,2);
z = X(:,3);

% Pressure 
nOfPoints = numel(x);
dp = zeros(nOfPoints,3);
dp(:,1) = 1-2*x;
dp(:,2) = 1-2*y;
dp(:,3) = 1-2*z;