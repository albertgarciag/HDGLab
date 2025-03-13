function du2 = poisson_3D_Ex1_analyticalD2(X,problemParams,iMat)

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

z1 = alpha*a*coA - beta*c*siC;
z2 = alpha*b*coA - beta*d*siC;
z3 = alpha*e*coA - beta*f*siC;

% Solution
u = exp(alpha*siA + beta*coC);
ux = u.*z1;
uy = u.*z2;
uz = u.*z3;
uxx = ux.*z1 - u.*(alpha*a*a*siA + beta*c*c*coC);
uxy = uy.*z1 - u.*(alpha*a*b*siA + beta*c*d*coC);
uxz = uz.*z1 - u.*(alpha*a*c*siA + beta*c*e*coC);
uyy = uy.*z2 - u.*(alpha*b*b*siA + beta*d*d*coC);
uyz = uz.*z2 - u.*(alpha*b*c*siA + beta*d*e*coC);
uzz = uz.*z3 - u.*(alpha*e*e*siA + beta*f*f*coC);

du2 = [uxx, uxy, uxz, uxy, uyy, uyz, uxz, uyz, uzz];
