function t = stokes_slip(X, n, t, problemParams, iMat, nsd, example)

if nsd==2
    switch example
        case 1
            t = stokes_2D_Ex1_slip(X, n, t, problemParams, iMat);
        case 2
            t = stokes_2D_Ex2_slip(X, n, t, problemParams, iMat);
        case 4
            t = stokes_2D_Ex4_slip(X, n, t, problemParams, iMat);
    end
elseif nsd==3
    switch example
        case 1
            t = stokes_3D_Ex1_slip(X, n, t, problemParams, iMat);
        case 2
            t = stokes_3D_Ex2_slip(X, n, t, problemParams, iMat);
    end
end