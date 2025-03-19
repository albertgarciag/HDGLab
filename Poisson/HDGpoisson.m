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
setpathPoisson

%% User defined parameters
meshFile = 'squareH4P1';

hdg.tau = 1;
hdg.problem = 'Poisson';
problemParams.conductivity = 1;
problemParams.charLength = 1;
problemParams.example = 1;

outputPath = 'resPoisson';
computeError = 1;  

%% Computation
load(meshFile);
problemParams.nOfMat = max(mesh.matElem);

disp('HDG preprocess...')
[refElem, refFace] = getRefData(mesh.pElem, mesh.nsd, mesh.optionNodes);
ctt = hdg_Poisson_Constants();
[mesh,hdg] = hdgPreprocess(mesh,refElem,refFace,hdg,ctt);

disp('HDG solution...')
[uHat,local] = hdg_Poisson_GlobalSystem(mesh,refElem,refFace,hdg,ctt,problemParams);
[u,q] = hdg_Poisson_LocalProblem(mesh,refElem,refFace,hdg,uHat,local,ctt.nOfComponents);
uStar = hdg_Poisson_LocalPostprocess(mesh,refElem,u,q,problemParams,ctt.nOfComponents);
q = hdgReshape(q,ctt.nOfComponents*mesh.nsd);
gradU = hdg_Poisson_gradient(mesh, q, problemParams);

disp('Computing error...')
if computeError==1
    uErr = computeErrorL2norm(mesh,refElem,u,problemParams,ctt.nOfComponents,@poisson_analytical,uStar);
    qErr = computeErrorL2norm(mesh,refElem,gradU,problemParams,ctt.nOfComponents*mesh.nsd,@poisson_analyticalD,[]);
else    
    uErr=[];
    qErr=[];
end

solutionFile = sprintf('%s/%s_ex%d', outputPath, meshFile, problemParams.example);
save(solutionFile,'mesh','refElem','refFace','hdg','ctt','u','uHat','uStar','q','uErr','qErr');

%% Postprocess
hdg_Poisson_Postprocess