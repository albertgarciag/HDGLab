function a = convdiff_2D_Ex2_convection(X, problemParams, iMat)

% 30º constante
a = ones(size(X,1),2).*[cos(pi/6) sin(pi/6)];

% (x+1, y+1)
% a = X + [1 1];