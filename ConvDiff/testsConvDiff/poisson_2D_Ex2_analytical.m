function u = poisson_2D_Ex2_analytical(X,problemParams,iMat)

% Parameters
lambda = 10;
f = 6;

% Points
x = X(:,1);
y = X(:,2);

% Solution
u = 4*y.^2 - 4*lambda^2*y.*exp(-lambda*y).*cos(f*pi*x) + lambda*exp(-2*lambda*y);


