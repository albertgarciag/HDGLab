function du2 = poisson_2D_Ex2_analyticalD2(X,problemParams,iMat)

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
uxx = 4*pi^2*f^2*lambda^2*y.*z1.*c;
uxy = -4*pi*f*lambda^2*sin(pi*f*x).*z1.*(lambda*y - 1);
uyy = 4*lambda^3*exp(-2*lambda*y) + 8*lambda^3*z1.*c - 4*lambda^4*y.*z1.*c + 8;
du2 = [uxx, uxy, uxy, uyy];


