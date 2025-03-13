function u = convdiff_convection(X,problemParams,iMat,nsd,example)

if nsd==2
    switch example
        case 1
            u = convdiff_2D_Ex1_convection(X,problemParams,iMat);
        case 2
            u = convdiff_2D_Ex2_convection(X,problemParams,iMat);
    end
elseif nsd==3
    switch example
        case 1
            disp('Example not implemented')
            u = convdiff_3D_Ex1_convection(X,problemParams,iMat);
        case 2
            disp('Example not implemented')
            u = convdiff_3D_Ex2_convection(X,problemParams,iMat);
    end
end