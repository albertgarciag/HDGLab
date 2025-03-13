function u = poisson_Dirichlet(X,problemParams,iMat,nsd,example)

if nsd==2
    switch example
        case 1
            u = poisson_2D_Ex1_Dirichlet(X,problemParams,iMat);
        case 2
            u = poisson_2D_Ex2_Dirichlet(X,problemParams,iMat);
    end
elseif nsd==3
    switch example
        case 1
            u = poisson_3D_Ex1_Dirichlet(X,problemParams,iMat);
        case 2
            u = poisson_3D_Ex2_Dirichlet(X,problemParams,iMat);
    end
end