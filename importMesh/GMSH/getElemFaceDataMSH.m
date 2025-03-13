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
%
% [elemID, faceID, lowerEntities, nOfElemNodes, nOfFaceNodes] = getElemFaceDataMSH(pDegree, nsd)
%
% Inputs:
%
% pDegree:        polynomial degree of approximation
% nsd:            number of spatial dimensions
%
% Outputs:
%
% elemID:         GMSH ID of the element type
% faceID:         GMSH ID of the face type
% lowerEntities:  GMSH IDs of the geometrical entities of lower dimensions
%                 (points in 2D and points&lines in 3D)
% nOfElemNodes:   number of element nodes
% nOfFaceNodes:   number of face nodes
%

function [elemID, faceID, lowerEntities, nOfElemNodes, nOfFaceNodes] = getElemFaceDataMSH(pDegree, nsd)

if nsd == 2 %%%%%%%%%%%%%%%%%%%%% TRI
   lowerEntities = [15]; % points
   switch pDegree
       case 1
           faceID = 1;
           nOfFaceNodes = 2;
           elemID = 2;
           nOfElemNodes = 3;
       case 2
           faceID = 8;
           nOfFaceNodes = 3;
           elemID = 9;
           nOfElemNodes = 6;
       case 3
           faceID = 26;
           nOfFaceNodes = 4;
           elemID = 21;
           nOfElemNodes = 10;
       case 4
           faceID = 27;
           nOfFaceNodes = 5;
           elemID = 23;
           nOfElemNodes = 15;
       case 5
           faceID = 28;
           nOfFaceNodes = 6;
           elemID = 25;
           nOfElemNodes = 21;
   end
else %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TET
   lowerEntities = [15,1,8,26,27,28]; % points & lines
   switch pDegree
       case 1
           faceID = 2;
           nOfFaceNodes = 3;
           elemID = 4;
           nOfElemNodes = 4;
       case 2
           faceID = 9;
           nOfFaceNodes = 6;
           elemID = 11;
           nOfElemNodes = 10;
       case 3
           faceID = 21;
           nOfFaceNodes = 10;
           elemID = 29;
           nOfElemNodes = 20;
       case 4
           faceID = 23;
           nOfFaceNodes = 15;
           elemID = 30;
           nOfElemNodes = 35;
       case 5
           faceID = 25;
           nOfFaceNodes = 21;
           elemID = 31;
           nOfElemNodes = 56;
   end
end