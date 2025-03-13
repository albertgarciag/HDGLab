function t= convdiff_2D_Ex2_Neumann(X, n, problemParams, iMat)

% NOTE that the normal is taken as an input so
% this is only valid for polygonal boundaries or NEFEM (exact normal)
% otherwise t should be defined analytically

t = zeros(size(X,1),1);