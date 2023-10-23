function res = minus(obj1,obj2)
% minus - overloads the - operator for stl objects
%
% Syntax:
%    res = minus(obj1,obj2)
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
%    eq = x(1) - x(2) < 5
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

    % different cases depending on the the type of inputs 
    if isa(obj1,'stl') && ~obj1.logic
        if isa(obj2,'stl') && ~obj2.logic
            
            res = obj1;
            
            res.type = '-';
            res.lhs = obj1;
            res.rhs = obj2;
            res.variables = unique([obj1.variables;obj2.variables]);
            
        elseif isnumeric(obj2) && isscalar(obj2)
            
            res = obj1;
            
            res.type = '-';
            res.lhs = obj1;
            res.rhs = obj2;
            
        else       
            throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
        end
        
    elseif isa(obj2,'stl') && ~obj2.logic && isnumeric(obj1) && isscalar(obj1)
        
        res = obj2;
        
        res.type = '-';
        res.lhs = obj1;
        res.rhs = obj2;
        
    else
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
end

% ------------------------------ END OF CODE ------------------------------
