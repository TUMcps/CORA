function [p] = randPointExtreme(obj)
% randPointExtreme - generates a random extreme point of a zonotope bundle
%
% Syntax:  
%    [p] = randPointExtreme(obj)
%
% Inputs:
%    obj - zonotope bundle object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% obtain vertices
Vmat = vertices(obj);

% obtain number of vertices
nrOfVertices = length(Vmat(1,:));

% random vertex
randInd = ceil(rand*nrOfVertices);

% random point
p = Vmat(:,randInd);

%------------- END OF CODE --------------