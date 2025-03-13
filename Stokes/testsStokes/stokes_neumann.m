function t = stokes_neumann(X, n, problemParams, iMat, nsd, example)

if nsd==2
    switch example
        case 1
            t = stokes_2D_Ex1_neumann(X, n, problemParams, iMat);
        case 2
            t = stokes_2D_Ex2_neumann(X, n, problemParams, iMat);
        case 4
            t = stokes_2D_Ex4_neumann(X, n, problemParams, iMat);
    end
elseif nsd==3
    switch example
        case 1
            t = stokes_3D_Ex1_neumann(X, n, problemParams, iMat);
        case 2
            t = stokes_3D_Ex2_neumann(X, n, problemParams, iMat);
    end
end