function h = stokes_3D_Ex1_slip(X, n, t, problemParams,iMat)
% Synthetic solution in a unit cube

[ngp,nsd] = size(X);

v  = stokes_3D_Ex1_analyticalVelocity(X,problemParams,iMat);
dv = stokes_3D_Ex1_analyticalVelocityD(X,problemParams,iMat);
p  = stokes_3D_Ex1_analyticalPressure(X,problemParams,iMat);
nu = problemParams.viscosity(iMat);

D = [n, problemParams.betaSlip*t];
E = [problemParams.alphaSlip*n, t];

u1x = dv(:,1);
u1y = dv(:,2);
u1z = dv(:,3);
u2x = dv(:,4);
u2y = dv(:,5);
u2z = dv(:,6);
u3x = dv(:,7);
u3y = dv(:,8);
u3z = dv(:,9);

gradN = [u1x.*n(:,1) + u1y.*n(:,2) + u1z.*n(:,3), ...
         u2x.*n(:,1) + u2y.*n(:,2) + u2z.*n(:,3), ...
         u3x.*n(:,1) + u3y.*n(:,2) + u3z.*n(:,3)];
                                            
h = zeros(ngp,nsd);
for iNsd = 1:nsd
    for jNsd = 1:nsd
        h(:,iNsd) = h(:,iNsd) + v(:,jNsd).*D(:,nsd*(iNsd-1)+jNsd) + (nu*gradN(:,jNsd) - p.*n(:,jNsd)).*E(:,nsd*(iNsd-1)+jNsd);
    end
end