function u = stokes_analyticalVelocity(X, problemParams, iMat, nsd, example)

if nsd==2
    switch example
        case 1
            u = stokes_2D_Ex1_analyticalVelocity(X, problemParams, iMat);
        case 2
            u = stokes_2D_Ex2_analyticalVelocity(X, problemParams, iMat);
        case 3
            u = stokes_2D_Ex3_analyticalVelocity(X, problemParams, iMat);
        case 4
            u = stokes_2D_Ex4_analyticalVelocity(X, problemParams, iMat);
    end
elseif nsd==3
    switch example
        case 1
            u = stokes_3D_Ex1_analyticalVelocity(X, problemParams, iMat);
        case 2
            u = stokes_3D_Ex2_analyticalVelocity(X, problemParams, iMat);
    end
end