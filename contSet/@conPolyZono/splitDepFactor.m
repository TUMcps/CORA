function cPZsplit = splitDepFactor(obj,factor)
% splitDepFactor - Splits one factor of a conPolyZono object
%
% Syntax:  
%    cPZsplit = splitDepFactor(obj,factor)
%
% Inputs:
%    obj - conPolyZono object
%    factor - index of the factor whos domain is splitted
%
% Outputs:
%    cPZsplit - cell array of split conPolyZono objects

% Example: 
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split, splitLongestGen

% Author:       Niklas Kochdumper
% Written:      07-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % construct domains after split
    dom1 = interval(-1,0);
    dom2 = interval(0,1);
    
    % compute the splitted conPolyZonotope objects
    cPZsplit{1} = getSubset(obj,factor,dom1);
    cPZsplit{2} = getSubset(obj,factor,dom2);
    
end
    
%------------- END OF CODE --------------