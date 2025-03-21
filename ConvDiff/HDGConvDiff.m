%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2D/3D HDG code for the Poisson equation using isoparametric elements   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ruben Sevilla, Matteo Giacomini and Antonio Huerta    
% Contact:      r.sevilla@swansea.ac.uk                                   
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
setpathConvDiff

%% User defined parameters
meshFile = 'squareH3P4';
beta = 10.^linspace(-2,2.8,13);

hdg.tau_d = 1;
hdg.beta = 40;
hdg.problem = 'ConvDiff';
problemParams.conductivity = 1;
problemParams.charLength = 1;
problemParams.example = 2;
% 1 if NBC evaluates the total flux, 0 if only q
problemParams.totalFluxNeumann = 1;

outputPath = 'resConvDiff';
computeError = 0;  

%% Computation
load(meshFile);
problemParams.nOfMat = max(mesh.matElem);
mesh.extFaces(mesh.extFaces(:,4)==1,3)=1;
mesh.extFaces(mesh.extFaces(:,4)==2,3)=2;
mesh.extFaces(mesh.extFaces(:,4)==3,3)=2;
mesh.extFaces(mesh.extFaces(:,4)==4,3)=1;

%computation of tau_d and tau_a
a = convdiff_convection(mesh.X,problemParams,[],mesh.nsd,problemParams.example);
hdg.tau_a = hdg.beta * max(sqrt(sum(a.^2,2)));
hdg.tau = hdg.tau_d + hdg.tau_a;

disp('HDG preprocess...')
[refElem, refFace] = getRefData(mesh.pElem, mesh.nsd, mesh.optionNodes);
ctt = hdg_ConvDiff_Constants();
[mesh,hdg] = hdgPreprocess(mesh,refElem,refFace,hdg,ctt);

disp('HDG solution...')
[uHat,local] = hdg_ConvDiff_GlobalSystem(mesh,refElem,refFace,hdg,ctt,problemParams);
[u,q] = hdg_ConvDiff_LocalProblem(mesh,refElem,refFace,hdg,uHat,local,ctt.nOfComponents);
uStar = hdg_Poisson_LocalPostprocess(mesh,refElem,u,q,problemParams,ctt.nOfComponents);
q = hdgReshape(q,ctt.nOfComponents*mesh.nsd);
gradU = hdg_Poisson_gradient(mesh, q, problemParams);

disp('Computing error...')
if computeError==1
    uErr = computeErrorL2norm(mesh,refElem,u,problemParams,ctt.nOfComponents,@convdiff_analytical,uStar);
    qErr = computeErrorL2norm(mesh,refElem,gradU,problemParams,ctt.nOfComponents*mesh.nsd,@convdiff_analyticalD,[]);
else    
    uErr=[];
    qErr=[];
end

solutionFile = sprintf('%s/%s_ex%d_beta%d', outputPath, meshFile, problemParams.example,floor(hdg.beta*100));
save(solutionFile,'mesh','refElem','refFace','hdg','ctt','u','uHat','uStar','q','uErr','qErr');

%% Postprocess
hdg_Poisson_Postprocess
T = zeros(mesh.nOfElements, 3);
for iElem = 1:mesh.nOfElements
    T(iElem,:) = mesh.indexT(iElem,1):(mesh.indexT(iElem,1)+2);
end
figure(iFig+1)
trisurf(T, mesh.X(:,1), mesh.X(:,2), u)