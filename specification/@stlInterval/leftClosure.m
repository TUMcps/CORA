function int = leftClosure(obj)
% leftClosure - makes the interval left-closed
%
% Syntax:
%    int = leftClosure(obj)
%
% Inputs:
%    obj - stlInterval object
%
% Outputs:
%    int - left-closed stlInterval object
%
% Example:
%    int = stlInterval(1,2,false,true);
%    int = leftClosure(int);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: rightClosure

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isemptyobject(obj) || obj.lower == -inf
    int = obj;
else
    int = stlInterval(obj.lower,obj.upper,true,obj.rightClosed);
end

% ------------------------------ END OF CODE ------------------------------
