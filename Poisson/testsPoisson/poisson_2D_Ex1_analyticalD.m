function du = poisson_2D_Ex1_analyticalD(X,problemParams,iMat)

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

% Tmp
sA = sin(a*x+b*y);
sC = sin(c*x+d*y);
cA = cos(a*x+b*y);
cC = cos(c*x+d*y);

% Solution
u = exp(alpha*sA + beta*cC);

du = [u.*(alpha*a*cA - beta*c*sC), u.*(alpha*b*cA - beta*d*sC)];