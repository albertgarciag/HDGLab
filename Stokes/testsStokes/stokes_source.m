function s = stokes_source(X, problemParams, iMat, nsd, example)

if nsd==2
    switch example
        case 1
            s = stokes_2D_Ex1_source(X, problemParams, iMat);
        case 2
            s = stokes_2D_Ex2_source(X, problemParams, iMat);
        case 3
            s = stokes_2D_Ex3_source(X, problemParams, iMat);
        case 4
            s = stokes_2D_Ex4_source(X, problemParams, iMat);
    end
elseif nsd==3
    switch example
        case 1
            s = stokes_3D_Ex1_source(X, problemParams, iMat);
        case 2
            s = stokes_3D_Ex2_source(X, problemParams, iMat);
    end
end