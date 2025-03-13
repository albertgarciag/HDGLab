function p = stokes_analyticalPressureD(X, problemParams, iMat, nsd, example)

if nsd==2
    switch example
        case 1
            p = stokes_2D_Ex1_analyticalPressureD(X, problemParams, iMat);
        case 2
            p = stokes_2D_Ex2_analyticalPressureD(X, problemParams, iMat);
        case 4
            p = stokes_2D_Ex4_analyticalPressureD(X, problemParams, iMat);
    end
elseif nsd==3
    switch example
        case 1
            p = stokes_3D_Ex1_analyticalPressureD(X, problemParams, iMat);
        case 2
            p = stokes_3D_Ex2_analyticalPressureD(X, problemParams, iMat);
    end
end