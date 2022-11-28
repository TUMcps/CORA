function res = mtimes(obj1,obj2)
% mtimes - overloads the * operator for stl objects
%
% Syntax:  
%    res = mtimes(obj1,obj2)
%
% Inputs:
%    obj1 - first logic formula (class stl)
%    obj2 - second logic formula (class stl)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    eq = 2*x(2) < 5
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % different cases depending on the the type of inputs 
    if isa(obj1,'stl') && ~obj1.logic && isnumeric(obj2) && isscalar(obj2)
            
        res = obj1;
        
        res.type = '*';
        res.lhs = obj2;
        res.rhs = obj1;
        
    elseif isa(obj2,'stl') && ~obj2.logic && isnumeric(obj1) && isscalar(obj1)
        
        res = obj2;
        
        res.type = '*';
        res.lhs = obj1;
        res.rhs = obj2;
        
    else
        throw(CORAerror('CORA:notSupported',...
                      'This operation is not supported for stl objects!'));
    end
end

%------------- END OF CODE --------------