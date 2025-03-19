function du2 = convdiff_2D_Ex3_analyticalD2(X,problemParams,iMat)

% Parameters
% alpha = pi/30;
% beta = 3;
% 
% a = 7*pi/2;
% b = -2*pi;
% c = 3*pi;
% d = 2*pi;

alpha = 0.1;
beta = 0.3;

a = 5.1;
b = -6.2;
c = 4.3;
d = 3.4;

% Points
x = X(:,1);
y = X(:,2);
nOfPoints = numel(x);

% Tmp
sA = sin(a*x+b*y);
sC = sin(c*x+d*y);
cA = cos(a*x+b*y);
cC = cos(c*x+d*y);

z1 = alpha*a*cA - beta*c*sC;
z2 = alpha*b*cA - beta*d*sC;

% Solution
u = exp(alpha*sA + beta*cC);
du = [u.*z1, u.*z2];

du2 = zeros(nOfPoints,4);

du2(:,1) = du(:,1).*z1 - u.*( alpha*a^2*sA + beta*c^2*cC );
du2(:,2) = du(:,2).*z1 - u.*( alpha*a*b*sA + beta*c*d*cC );
du2(:,3) = du2(:,2);
du2(:,4) = du(:,2).*z2 - u.*( alpha*b^2*sA + beta*d^2*cC );
