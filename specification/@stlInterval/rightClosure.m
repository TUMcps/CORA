function int = rightClosure(obj)
% rightClosure - makes the interval right-closed
%
% Syntax:
%    int = rightClosure(obj)
%
% Inputs:
%    obj - stlInterval object
%
% Outputs:
%    int - right-closed stlInterval object
%
% Example:
%    int = stlInterval(1,2,true,false);
%    int = rightClosure(int);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: leftClosure

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isemptyobject(obj) || obj.upper == inf
    int = obj;
else
    int = stlInterval(obj.lower,obj.upper,obj.leftClosed,true);
end

% ------------------------------ END OF CODE ------------------------------
