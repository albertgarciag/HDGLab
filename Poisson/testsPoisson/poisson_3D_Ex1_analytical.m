function u = poisson_3D_Ex1_analytical(X,problemParams,iMat)

% Parameters
alpha = 0.1;
beta = 0.3;

a = 5.1;
b = -6.2;
c = 4.3;
d = 3.4;
e = 1.8;
f = 1.7;

% Points
x = X(:,1);
y = X(:,2);
z = X(:,3);

% Tmp
si = sin(a*x+b*y+e*z);
co = cos(c*x+d*y+f*z);

% Solution
u = exp(alpha*si + beta*co);