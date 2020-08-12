classdef matPolytope
% matPolytope class 
%
% Syntax:  
%    obj = matPolytope(vert)
%
% Inputs:
%    vert - cell-array storing the vertices
%
% Outputs:
%    obj - generated object
%
% Example:
%    V{1} = [1 2; 0 1];
%    V{2} = [1 3; -1 2];
%
%    mp = matPolytope(V);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix

% Author:       Matthias Althoff
% Written:      21-June-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    dim = 1;
    verts = 0; % number of vertices
    vertex = [];
end
    
methods
    
    % class constructor
    function obj = matPolytope(input)
        if nargin==1
            if isa(input, 'polytope')
                %get vertices from polytope class
                V=extreme(input);
                %rewrite vertices as matrices
                for i=1:length(V(:,1))
                    matrixVertex{i}=vec2mat(V(i,:));
                end
            else
                matrixVertex=input;
            end
            obj.dim = length(matrixVertex{1});
            obj.verts = length(matrixVertex);
            obj.vertex = matrixVertex;
        end
    end
         
    %methods in seperate files     
    matP = plus(summand1,summand2)
    matP = simplePlus(summand1,summand2)
    matP = mtimes(factor1,factor2)
    matP = mpower(matP,exponent)
    matP = expmInd(matP,maxOrder)
    [eP,eI] = expmIndMixed(matP,intermediateOrder,maxOrder)
    intMat = intervalMatrix(matP)
    matZ = matZonotope(matP)
    matP = polytope(matP)
    dist = expmDist(matP,exactMat,maxOrder)
    matV = vertices(matP)
    
    %display functions
    plot(varargin)
    plotComb(varargin)
    display(obj)
end
end

%------------- END OF CODE --------------