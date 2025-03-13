function t= poisson_3D_Ex2_Neumann(X, n, problemParams, iMat)

TOL = 1e-6;
xLeft = -5;

n = size(X,1);

t = zeros(n,1);
t(X(:,1)<xLeft+TOL) = -1;