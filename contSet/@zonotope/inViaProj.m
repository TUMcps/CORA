function result = inViaProj(Z1,Z2)
% inViaProj - determines if zonotope Z1 is enclosed by zonotope Z2 for all
% possible 2D projections
%
% Syntax:  
%    result = inViaProj(Z1,Z2)
%
% Inputs:
%    Z1 - 1st zonotope object
%    Z2 - 2nd zonotope object
%
% Outputs:
%    result - 1/0 if zonotope is in, or not
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
% Written:      06-April-2017 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% dimension
dim = length(center(Z1));

% combinations of projections
comb = combinator(dim,2,'c');

% init result
result = 1;

% loop through projections
for i = 1:length(comb(:,1))
    % perform projections
    Z1proj = project(Z1,comb(i,:));
    Z2proj = project(Z2,comb(i,:));
    
    % check enclosure
    if ~in(Z1proj, Z2proj)
        result = 0;
        break
    end
end


%------------- END OF CODE --------------