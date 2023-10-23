function res = implies(obj1,obj2)
% implies - implication operator for logic formulas
%
% Syntax:
%    res = implies(obj1,obj2)
%
% Inputs:
%    obj1 - first logic formula (class stl)
%    obj2 - second logic formula (class stl)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = implies(x(1) < 5,x(2) < 3)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check input arguments
    if ~isa(obj1,'stl') || ~isa(obj2,'stl') || ~obj1.logic || ~obj2.logic
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
    
    % construct resulting stl object
    res = obj1;
    
    res.type = '|';
    res.lhs = ~obj1;
    res.rhs = obj2;
    res.temporal = obj1.temporal | obj2.temporal;
    res.logic = obj1.logic & obj2.logic;
    res.from = [];
    res.to = [];
    res.variables = unique([obj1.variables;obj2.variables]);
end

% ------------------------------ END OF CODE ------------------------------
