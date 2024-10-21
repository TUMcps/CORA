function int = toLeft(obj)
% toLeft - get the interval to the negative side of the current one
%
% Syntax:
%    int = toLeft(obj)
%
% Inputs:
%    obj - stlInterval object
%
% Outputs:
%    int - stlInterval to the left of obj
%
% Example:
%    int = stlInterval(1,2,false,true);
%    int = toLeft(int);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: toRight

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
if isemptyobject(obj)
    throw(CORAerror('CORA:emptySet'));
end
% STL intervals are always bounded at 0
int = stlInterval(0,obj.lower,true,~obj.leftClosed);

% ------------------------------ END OF CODE ------------------------------
