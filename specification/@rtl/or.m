function res = or(obj1,obj2)
% or - overloads the | operator representing a logic "or" for rtl objects
%
% Syntax:
%    res = or(obj1,obj2)
%
% Inputs:
%    obj1 - first logic formula (class rtl)
%    obj2 - second logic formula (class rtl)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    eq = rtl(x(1) < 5) | rtl(x(2) < 3)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% quick exit
if isempty(obj1)
    res = obj2; return;
elseif isempty(obj2)
    res = obj1; return;
end

% check input arguments
inputArgsCheck({{obj1,'att','rtl'},{obj2,'att','rtl'}});

% construct resulting rtl object
res = obj1;

res.type = '|';
res.lhs = obj1;
res.rhs = obj2;
res.time = [];

% ------------------------------ END OF CODE ------------------------------
