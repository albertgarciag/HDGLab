function p = stokes_3D_Ex1_analyticalPressure(X,problemParams,iMat) 
% Synthetic solution in a unit cube

% Points
x = X(:,1);
y = X(:,2);
z = X(:,3);

% Pressure 
p = x.*(1-x) + y.*(1-y) + z.*(1-z);