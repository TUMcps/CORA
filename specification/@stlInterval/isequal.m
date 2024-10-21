function res = isequal(obj,other,tol)
% isequal - checks if two STL intervals are equal
%
% Syntax:
%    res = isequal(obj,other)
%    res = isequal(obj,other,tol)
%
% Inputs:
%    obj - stlInterval object
%    other - stlInterval object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    I1 = stlInterval(1,2,true,true);
%    I2 = stlInterval(1,2,false,true);
%    res = isequal(I1,I2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
narginchk(2,3);
% for speed:
if nargin < 3
    tol = eps;
end

inputArgsCheck({{obj,'att','stlInterval'};
                {other,'att','stlInterval'};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

if isemptyobject(obj) || isemptyobject(other)
    res = isemptyobject(obj) && isemptyobject(other);
else
    res = withinTol(obj.lower,other.lower,tol) && withinTol(obj.upper,other.upper,tol) && ...
        obj.leftClosed == other.leftClosed && obj.rightClosed == other.rightClosed;
end

% ------------------------------ END OF CODE ------------------------------
