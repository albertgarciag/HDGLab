function u = stokes_analyticalVelocityD2(X, problemParams, iMat, nsd, example)

if nsd==2
    switch example
        case 1
            u = stokes_2D_Ex1_analyticalVelocityD2(X, problemParams, iMat);
        case 2
            u = stokes_2D_Ex2_analyticalVelocityD2(X, problemParams, iMat);
        case 4
            u = stokes_2D_Ex4_analyticalVelocityD2(X, problemParams, iMat);
    end
elseif nsd==3
    switch example
        case 1
            u = stokes_3D_Ex1_analyticalVelocityD2(X, problemParams, iMat);
        case 2
            u = stokes_3D_Ex2_analyticalVelocityD2(X, problemParams, iMat);
    end
end