function res = and(obj1,obj2)
% and - overloads the & operator representing a logic "and" for rtl objects
%
% Syntax:
%    res = and(obj1,obj2)
%
% Inputs:
%    obj1 - first logic formula (class rtl)
%    obj2 - second logic formula (class rtl)
%
% Outputs:
%    res - resulting stl formula (class rtl)
%
% Example: 
%    x = stl('x',2)
%    eq = rtl(x(1) < 5) & rtl(x(2) < 3)
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

% check input arguments
if isempty(obj1)
    res = obj2; return;
elseif isempty(obj2)
    res = obj1; return;
end

% check input arguments
inputArgsCheck({{obj1,'att','rtl'},{obj2,'att','rtl'}});

% construct resulting rtl object
res = obj1;

res.type = '&';
res.lhs = obj1;
res.rhs = obj2;
res.time = [];

% ------------------------------ END OF CODE ------------------------------
