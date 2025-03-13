function h = stokes_2D_Ex4_slip(X, n, t, problemParams,iMat)

[ngp,nsd] = size(X);

v  = stokes_2D_Ex4_analyticalVelocity(X,problemParams,iMat);
dv = stokes_2D_Ex4_analyticalVelocityD(X,problemParams,iMat);
p  = stokes_2D_Ex4_analyticalPressure(X,problemParams,iMat);
nu = problemParams.viscosity(iMat);

D = [n, problemParams.betaSlip*t];
E = [problemParams.alphaSlip*n, t];

u1x = dv(:,1);
u1y = dv(:,2);
u2x = dv(:,3);
u2y = dv(:,4);

gradN = [u1x.*n(:,1) + u1y.*n(:,2), ...
         u2x.*n(:,1) + u2y.*n(:,2)];
                                            
h = zeros(ngp,nsd);
for iNsd = 1:nsd
    for jNsd = 1:nsd
        h(:,iNsd) = h(:,iNsd) + v(:,jNsd).*D(:,nsd*(iNsd-1)+jNsd) + (nu*gradN(:,jNsd) - p.*n(:,jNsd)).*E(:,nsd*(iNsd-1)+jNsd);
    end
end