function res = isequal(obj1,obj2)
% isequal - check if two STL objects are identical
%
% Syntax:
%    res = isequal(obj1,obj2)
%
% Inputs:
%    obj1 - logic formula (class stl)
%    obj2 - logic formula (class stl)
%
% Outputs:
%    res - identical (res = true) or not (res = false)
%
% Example: 
%    x = stl('x',2);
%    eq1 = x(1) + x(2) < 5 & x(2) > 2;
%    eq2 = x(2) > 2 & x(2) + x(1) < 5;
%    isequal(eq1,eq2)
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
    
    res = false;

    if isa(obj1,'stl') && isa(obj2,'stl') && strcmp(obj1.type,obj2.type)

        if ismember(obj1.type,{'true','false'})

            res = true;

        elseif strcmp(obj1.type,'variable')

            res = strcmp(obj1.var,obj2.var);

        elseif ismember(obj1.type,{'&','|','+','*'})

            if isequal(obj1.lhs,obj2.lhs)
                res = isequal(obj1.rhs,obj2.rhs);
            else
                if isequal(obj1.lhs,obj2.rhs)
                    res = isequal(obj1.rhs,obj2.lhs);
                end
            end

        elseif ismember(obj1.type,{'-','<','<=','>','>='})

            if isequal(obj1.lhs,obj2.lhs) && isequal(obj1.rhs,obj2.rhs)
                res = true;
            end

        elseif ismember(obj1.type,{'finally','globally','next'})

            if isequal(obj1.from,obj2.from) && isequal(obj1.to,obj2.to)
                res = isequal(obj1.lhs,obj2.lhs);
            end

        elseif ismember(obj1.type,{'until','release'})

            if isequal(obj1.from,obj2.from) && isequal(obj1.to,obj2.to) ...
                && isequal(obj1.lhs,obj2.lhs) && isequal(obj1.rhs,obj2.rhs)
                res = true;
            end

        else
            
            res = isequal(obj1.lhs,obj2.lhs);           
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
