function a = convdiff_2D_Ex1_source(X, problemParams, iMat)

a = [-sin(pi*X(:,1)).*cos(pi*X(:,2)), cos(pi*X(:,1)).*sin(pi*X(:,2))];