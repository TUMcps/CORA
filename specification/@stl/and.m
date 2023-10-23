function res = and(obj1,obj2)
% and - overloads the & operator representing a logic "and" for stl objects
%
% Syntax:
%    res = and(obj1,obj2)
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
%    eq = x(1) < 5 & x(2) < 3
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
    if isempty(obj1)
        res = obj2; return;
    elseif isempty(obj2)
        res = obj1; return;
    end

    obj1 = aux_bool2stl(obj1); obj2 = aux_bool2stl(obj2);

    if ~isa(obj1,'stl') || ~isa(obj2,'stl') || ~obj1.logic || ~obj2.logic
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
    
    % construct resulting stl object
    res = obj1;
    
    res.type = '&';
    res.lhs = obj1;
    res.rhs = obj2;
    res.temporal = obj1.temporal | obj2.temporal;
    res.logic = obj1.logic & obj2.logic;
    res.from = [];
    res.to = [];
    res.id = [];
    res.variables = unique([obj1.variables;obj2.variables]);
end


% Auxiliary functions -----------------------------------------------------

function obj = aux_bool2stl(obj)
% convert a boolean value to a stl object

    if islogical(obj)
        if obj
            obj = stl('x',1);
            obj.type = 'true';
        else
            obj = stl('x',1);
            obj.type = 'false';
        end
        obj.variables = [];
        obj.logic = true;
    end
end

% ------------------------------ END OF CODE ------------------------------
