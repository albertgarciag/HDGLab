%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2D/3D HDG code for the Stokes equation using isoparametric elements   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ruben Sevilla, Matteo Giacomini and Antonio Huerta    
% Contact:      matteo.giacomini@upc.edu                                   
% ZCCE (SU)   - https://www.swansea.ac.uk/engineering/zcce/
% LaCaN (UPC) - https://www.lacan.upc.edu/                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc
close all
setpathStokes

%% User defined parameters
meshFile = 'squareH1P2';
% meshFile = 'cubeH3P5';
% meshFile = 'ringH4P5';
% meshFile = 'sphereExtP5';

hdg.tau = 3;
hdg.problem = 'Stokes';
problemParams.viscosity = 1;
problemParams.alphaSlip = 0;
problemParams.betaSlip = 0;
problemParams.charLength = 1;
problemParams.example = 4;
problemParams.axi = 1; 

outputPath = 'resStokes';
computeError = 1;  

%% Computation
load(meshFile);
problemParams.nOfMat = max(mesh.matElem);    

disp('HDG preprocess...')
[refElem, refFace] = getRefData(mesh.pElem, mesh.nsd, mesh.optionNodes);
ctt = hdg_Stokes_Constants(mesh.nsd);

% Change boundary conditions
% 1: all Dirichlet
% 2: replace Neumann by slip BC
% 3: sphere for external flow
% 4: axisymmetric Poiseuille (bottom horizontal axis of symmetry)
changeBC = 4; 
switch changeBC
    case 1
        mesh.extFaces(mesh.extFaces(:,3)==ctt.iBC_Neumann,3) = ctt.iBC_Dirichlet;
    case 2
        mesh.extFaces(mesh.extFaces(:,3)==ctt.iBC_Neumann,3) = ctt.iBC_Slip;
    case 3
        % Top
        mesh.extFaces(mesh.extFaces(:,4)==1,3) = ctt.iBC_Dirichlet;
        % Inlet
        mesh.extFaces(mesh.extFaces(:,4)==2,3) = ctt.iBC_Dirichlet;
        % Outlet
        mesh.extFaces(mesh.extFaces(:,4)==3,3) = ctt.iBC_Neumann;
        % Lateral 
        mesh.extFaces(mesh.extFaces(:,4)==4,3) = ctt.iBC_Dirichlet;
        % Sphere
        mesh.extFaces(mesh.extFaces(:,4)==5,3) = ctt.iBC_Dirichlet;
        mesh.extFaces(mesh.extFaces(:,4)==6,3) = ctt.iBC_Dirichlet;
        % Bottom
        mesh.extFaces(mesh.extFaces(:,4)==7,3) = ctt.iBC_Slip;
        % Lateral (sphere)
        mesh.extFaces(mesh.extFaces(:,4)==8,3) = ctt.iBC_Slip;
    case 4
        % Top
        mesh.extFaces(mesh.extFaces(:,4)==2,3) = ctt.iBC_Dirichlet;
        % Bottom
        mesh.extFaces(mesh.extFaces(:,4)==4,3) = ctt.iBC_Axi;
        % Left
        mesh.extFaces(mesh.extFaces(:,4)==1,3) = ctt.iBC_Dirichlet;
        % Right
        mesh.extFaces(mesh.extFaces(:,4)==3,3) = ctt.iBC_Neumann;
    otherwise
        disp('No change to default BC')
end

[mesh, hdg] = hdgPreprocess(mesh, refElem, refFace, hdg, ctt);

disp('HDG solution...')
[uHat, rho, local] = hdg_Stokes_GlobalSystem(mesh, refElem, refFace, hdg, ctt, problemParams);
[u, p, L] = hdg_Stokes_LocalProblem(mesh, refElem, refFace, hdg, uHat, rho, local, ctt.nOfComponents);
uStar = hdg_Stokes_LocalPostprocess(mesh, refElem, u, L, problemParams, ctt.nOfComponents);

u = hdgReshape(u,ctt.nOfComponents);
uStar = hdgReshape(uStar,ctt.nOfComponents);
L = hdgReshape(L,ctt.nOfComponents*mesh.nsd);
gradU = hdg_Stokes_gradient(mesh, L, problemParams);

if computeError==1
    disp('Computing error...')
    uErr = computeErrorL2norm(mesh,refElem,u,problemParams,ctt.nOfComponents,@stokes_analyticalVelocity,uStar);
    pErr = computeErrorL2norm(mesh,refElem,p,problemParams,1,@stokes_analyticalPressure,[]);
    LErr = computeErrorL2norm(mesh,refElem,gradU,problemParams,ctt.nOfComponents*mesh.nsd,@stokes_analyticalVelocityD,[]);
else    
    uErr=[];
    pErr=[];
    LErr=[];
end

solutionFile = sprintf('%s/%s_ex%d', outputPath, meshFile, problemParams.example);
save(solutionFile,'mesh','refElem','refFace','hdg','ctt','u','uHat','uStar','p','rho','L','uErr','pErr','LErr');

%% Postprocess
hdg_Stokes_Postprocess