function u = convdiff_2D_Ex2_Dirichlet(X,problemParams,iMat)

l1 = X(:,1)==0;
l2 = X(:,2)<=0.2;
u = zeros(size(X,1),1);
u(logical(l1.*l2)) = 1;