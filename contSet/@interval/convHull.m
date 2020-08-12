function res = convHull(Int1,Int2)
% convHull - computes an enclosure for the convex hull of two intervals
%
% Syntax:  
%    res = convHull(Int1,Int2)
%
% Inputs:
%    Int1 - first interval object
%    Int2 - second interval object
%
% Outputs:
%    res - interval enclosing the convex hull of Int1 and Int2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/convHull

% Author:        Niklas Kochdumper
% Written:       26-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

%------------- BEGIN CODE --------------

% determine the interval object
if ~isa(Int1,'interval')
    temp = Int1;
    Int1 = Int2;
    Int2 = temp;
end

% different cases depending on the class of the summand
if isa(Int2,'interval') || isnumeric(Int2)

    res = Int1 | Int2;

elseif isa(Int2,'zonotope') || isa(Int2,'conZonotope') || ...
       isa(Int2,'zonoBundle') || isa(Int2,'polyZonotope') || ...
       isa(Int2,'mptPolytope')

    res = convHull(Int2,Int1);

else
    
    % throw error for given arguments
    error(noops(Int1,Int2));
        
end

%------------- END OF CODE --------------