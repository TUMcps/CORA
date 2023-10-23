function cPZsplit = splitDepFactor(cPZ,factor)
% splitDepFactor - Splits one factor of a constrained polynomial zonotope
%
% Syntax:
%    cPZsplit = splitDepFactor(cPZ,factor)
%
% Inputs:
%    cPZ - conPolyZono object
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

% Authors:       Niklas Kochdumper
% Written:       07-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{cPZ,'att','conPolyZono'};
                {factor,'att','numeric',{'nonnan','integer',...
                    'nonnegative','scalar'}}});

% construct domains after split
dom1 = interval(-1,0);
dom2 = interval(0,1);

% compute the splitted conPolyZonotope objects
cPZsplit{1} = getSubset(cPZ,factor,dom1);
cPZsplit{2} = getSubset(cPZ,factor,dom2);

% ------------------------------ END OF CODE ------------------------------
