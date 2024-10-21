function int = toRight(obj)
% toRight - get the interval to the positive side of the current one
%
% Syntax:
%    int = toRight(obj)
%
% Inputs:
%    obj - stlInterval object
%
% Outputs:
%    int - stlInterval to the right of obj
%
% Example:
%    int = stlInterval(1,2,false,true);
%    int = toRight(int);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: toLeft

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isemptyobject(obj)
    throw(CORAerror('CORA:emptySet'));
end
int = stlInterval(obj.upper,inf,~obj.rightClosed,false);

% ------------------------------ END OF CODE ------------------------------
