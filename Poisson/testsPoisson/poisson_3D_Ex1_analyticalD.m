function du = poisson_3D_Ex1_analyticalD(X,problemParams,iMat)

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
siA = sin(a*x+b*y+e*z);
coC = cos(c*x+d*y+f*z);
coA = cos(a*x+b*y+e*z);
siC = sin(c*x+d*y+f*z);

% Solution
u = exp(alpha*siA + beta*coC);
ux = u.*(alpha*a*coA - beta*c*siC);
uy = u.*(alpha*b*coA - beta*d*siC);
uz = u.*(alpha*e*coA - beta*f*siC);
du = [ux, uy, uz];