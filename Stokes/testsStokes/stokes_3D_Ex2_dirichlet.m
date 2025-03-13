function v = stokes_3D_Ex2_dirichlet(X, problemParams, iMat) 
% External flow past a sphere

xInlet = -5;
zTop = 5;
yLat = 5;

tol = 1e-4;

% Initialise null velocity
v = zeros(size(X));

if all(X(:,1) < xInlet + tol) || all(X(:,2) > yLat - tol) || all(X(:,3) > zTop - tol)
    v = stokes_3D_Ex2_analyticalVelocity(X, problemParams, iMat);
end