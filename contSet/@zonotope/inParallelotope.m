function [res] = inParallelotope(Z1,Z2)
% inParallelotope - checks if a zonotope Z1 is in a parallelotope Z2, which
% is represented as a zonotope
%
% Syntax:  
%    [res] = inParallelotope(Z1,Z2)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - parallelotope in zonotope representation
%
% Outputs:
%    res - 1/0 depending if Z1 is enclosed in Z2 or not
%
% Example: 
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also: in

% Author:       Matthias Althoff
% Written:      23-September-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%init 
res = [];
dim = length(center(Z1));

%check if Z2 is a parallelotope
if isa(Z2,'zonotope')
    %obtain matrix of generators
    G = Z2.Z(:,2:end);
    %obtain number of generators
    nrOfGenerators = length(G(1,:));
    % Is Z2 a parallelotope?
    if (nrOfGenerators == dim) && (det(G) ~= 0)
        %obtain transformation matrix
        T = inv(G);
        
        %map Z1 and Z2
        Z1_trans = T*Z1;
        Z2_trans = T*Z2;
        
        %check enclosure
        res = interval(Z1_trans)<=interval(Z2_trans);
    end
end


%------------- END OF CODE --------------