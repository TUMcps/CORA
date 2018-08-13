function result = in(Z1,Z2)
% in - determines if zonotope Z1 is enclosed by zonotope Z2
%
% Syntax:  
%    result = in(Z1,Z2)
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
% Written:      07-May-2007 
% Last update:  06-April-2017
% Last revision:---

%------------- BEGIN CODE --------------

% generate halfspace representation if empty
if isempty(Z2.halfspace)
    Z2 = halfspace(Z2);
end

%simple test: Is the center of Z1 in Z2?
c = center(Z1);
inequality = (Z2.halfspace.H*c<=Z2.halfspace.K);

if ~all(inequality)
    result = 0;
else
    % perform halfspace projection
    halfspaceProj = Z2.halfspace.H*Z1 + (-Z2.halfspace.K);
    
    % obtain intervals of halfspace projection
    intHalfProj = interval(halfspaceProj);
    
    % supremum has to be less or equal to 0
    result = all(supremum(intHalfProj) <= 0);

end


%------------- END OF CODE --------------