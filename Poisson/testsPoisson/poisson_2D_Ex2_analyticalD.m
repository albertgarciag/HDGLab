function du = poisson_2D_Ex2_analyticalD(X,problemParams,iMat)

% Parameters
lambda = 10;
f = 6;

% Points
x = X(:,1);
y = X(:,2);

% Tmp
z1 = exp(-lambda*y);
c = cos(pi*f*x);

% Solution
ux = 4*pi*f*lambda^2*y.*sin(pi*f*x).*z1;
uy = 8*y - 2*lambda^2*exp(-2*lambda*y) - 4*lambda^2*z1.*c + 4*lambda^3*y.*z1.*c;
du = [ux, uy];


