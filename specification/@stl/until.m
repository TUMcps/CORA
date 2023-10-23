function res = until(obj1,obj2,time)
% until - until-operator for Signal Temporal Logic
%
% Syntax:
%    res = until(obj1,obj2,time)
%
% Inputs:
%    obj1 - first logic formula (class stl)
%    obj2 - second logic formula (class stl)
%    time - time interval (class interval)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = until(x(1) < 5, x(2) > 3,interval(0.1,0.2))
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
    obj1 = aux_bool2stl(obj1); obj2 = aux_bool2stl(obj2);

    if ~isa(obj1,'stl') || ~isa(obj2,'stl') || ~obj1.logic || ~obj2.logic
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
    
    if ~isa(time,'interval') || ~all(size(time) == [1,1])
        throw(CORAerror('CORA:wrongValue',...
                        'Wrong format for input argument "time"!'));
    end

    % construct resulting stl object
    res = obj1;
    
    res.type = 'until';
    res.lhs = obj1;
    res.rhs = obj2;
    res.id = [];
    res.temporal = true;
    res.from = infimum(time);
    res.to = supremum(time);
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
