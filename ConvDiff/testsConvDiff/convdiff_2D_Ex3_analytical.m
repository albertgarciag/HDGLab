function u = convdiff_2D_Ex3_analytical(X,problemParams,iMat)

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
s = sin(a*x+b*y);
c = cos(c*x+d*y);

% Solution
u = exp(alpha*s + beta*c);