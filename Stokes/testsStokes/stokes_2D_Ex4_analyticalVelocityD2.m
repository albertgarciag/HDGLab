function dv2 = stokes_2D_Ex4_analyticalVelocityD2(X, problemParams,iMat) 

% Data
R = 1;
V = 1;

% Points
nOfPoints = size(X,1);

% Second derivatives of the velocity
v1zz = zeros(nOfPoints,1);
v1zr = zeros(nOfPoints,1);
v1rz = zeros(nOfPoints,1);
v1rr = -2*V/R^2*ones(nOfPoints,1);
v2zz = zeros(nOfPoints,1);
v2zr = zeros(nOfPoints,1);
v2rz = zeros(nOfPoints,1);
v2rr = zeros(nOfPoints,1);

dv2 = [v1zz, v1zr, v1rz, v1rr, v2zz, v2zr, v2rz, v2rr];