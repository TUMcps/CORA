function res = combineNext(obj)
% combineNext - remove redundant next-opearators in an STL formula 
%
% Syntax:
%    res = combineNext(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    res - resulting logic formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = next(next(x(1) < 5 | x(2) < 3,0.1) | x(2) > 5,0.2);
%    eq_ = combineNext(eq)
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

    if ~obj.temporal
        res = obj;
    elseif strcmp(obj.type,'next')
        
        if strcmp(obj.lhs.type,'&') && obj.lhs.temporal
            res = combineNext(next(obj.lhs.lhs,obj.from)) & ...
                            combineNext(next(obj.lhs.rhs,obj.from));
        elseif strcmp(obj.lhs.type,'|') && obj.lhs.temporal
            res = combineNext(next(obj.lhs.lhs,obj.from)) | ...
                            combineNext(next(obj.lhs.rhs,obj.from));
        elseif strcmp(obj.lhs.type,'~')
            res = ~combineNext(next(obj.lhs.lhs,obj.from));
        elseif strcmp(obj.lhs.type,'next')
            res = combineNext(next(obj.lhs.lhs,obj.from + obj.lhs.from));
        else
            res = obj;
        end

    elseif strcmp(obj.type,'&')
        res = combineNext(obj.lhs) & combineNext(obj.rhs);
    elseif strcmp(obj.type,'|')
        res = combineNext(obj.lhs) | combineNext(obj.rhs);
    elseif strcmp(obj.type,'~')
        res = ~combineNext(obj.lhs);
    else
        res = obj;
    end  
end

% ------------------------------ END OF CODE ------------------------------
