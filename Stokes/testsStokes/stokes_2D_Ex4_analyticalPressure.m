function p = stokes_2D_Ex4_analyticalPressure(X,problemParams,iMat) 

% Data
R = 1;
V = 1;

% Points
z = X(:,1); 

% Pressure 
p = -4*V*z/R^2;