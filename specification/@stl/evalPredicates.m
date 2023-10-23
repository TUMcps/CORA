function res = evalPredicates(obj,timeInt,varargin)
% evalPredicates - replace predicates with there evaluation (true or false)
%
% Syntax:
%    res = evalPredicates(obj,timeInt)
%    res = evalPredicates(obj,timeInt,timePoint)
%
% Inputs:
%    obj - logic formula (class stl)
%    timeInt - matrix where the entry at index (i,j) speficies that the
%              predicate with id i at time step j is satisfied (value 1), 
%              violated (value 0), or unknown (value 2) for the time
%              interval reachable set
%    timePoint - same a timeInt but for the time point reachable set
%
% Outputs:
%    res - resulting logic formula (class stl)
%
% Example: 
%    x = stl('x',2);
%    eq = x(1) < 5 | globally(x(2) < 3,interval(1,2));
%    eq = assignIdentifiers(eq);
%    pred = [false false true; true 2 2];
%    res = evalPredicates(eq,pred)
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
    
    % initialization
    pred.timeInt = timeInt;

    if nargin > 2 && ~isempty(varargin{1})
        pred.timePoint = varargin{1};
    else
        pred.timePoint = 2*ones(size(timeInt,1),size(timeInt,2)+1);
    end

    % substitute all predicates by true/false using recursive function
    [res,mod] = aux_recursive(obj,1,pred);

    % eliminate true/false values from the resulting formula
    if mod
        res = eliminateTrueFalse(res);
    end
end


% Auxiliary functions -----------------------------------------------------

function [res,mod] = aux_recursive(obj,time,pred)
% recursive function to assign ids to all non-temporal formulas

    res = obj; mod = false;

    if ~obj.temporal

        if ~isempty(obj.id) && obj.id <= size(pred.timePoint,1)
            if pred.timePoint(obj.id,time) == 1
                res = stl(true); mod = true;
            elseif pred.timePoint(obj.id,time) == 0
                res = stl(false); mod = true;
            end
        end

    elseif strcmp(obj.type,'next')
        
        [res.lhs,mod] = aux_recursive(obj.lhs,obj.from+time,pred);

    elseif strcmp(obj.type,'globally') 

        if ~obj.lhs.temporal && ~isempty(obj.lhs.id) && ...
           obj.lhs.id <= size(pred.timePoint,1) && ...
           obj.lhs.id <= size(pred.timeInt,1) 
            
            timeInt = pred.timeInt(obj.lhs.id,time+obj.from:time+obj.to-1);
            timePoint = pred.timePoint(obj.lhs.id,time+obj.from:time+obj.to);

            if all(timeInt == 1)
                res = stl(true); mod = true;
            elseif any(timeInt == 0) || any(timePoint == 0)
                res = stl(false); mod = true;
            end
        end

    elseif strcmp(obj.type,'finally') 

        if ~obj.lhs.temporal && ~isempty(obj.lhs.id) && ...
           obj.lhs.id <= size(pred.timePoint,1) && ...
           obj.lhs.id <= size(pred.timeInt,1) 
            
            timeInt = pred.timeInt(obj.lhs.id,time+obj.from:time+obj.to-1);
            timePoint = pred.timePoint(obj.lhs.id,time+obj.from:time+obj.to);

            if any(timePoint == 1)
                res = stl(true); mod = true;
            elseif all(timeInt == 0)
                res = stl(false); mod = true;
            end
        end

    elseif strcmp(obj.type,'&')

        [res.lhs,mod1] = aux_recursive(obj.lhs,time,pred);

        if strcmp(res.lhs,'false')
            res = stl(false); mod = true;
        elseif strcmp(res.lhs,'true')
            [res,mod2] = aux_recursive(obj.rhs,time,pred);
            mod = mod1 | mod2;
        else
            [res.rhs,mod2] = aux_recursive(obj.rhs,time,pred);
            mod = mod1 | mod2;
        end

    elseif strcmp(obj.type,'|')

        [res.lhs,mod1] = aux_recursive(obj.lhs,time,pred);

        if strcmp(res.lhs,'true')
            res = stl(true); mod = true;
        elseif strcmp(res.lhs,'false')
            [res,mod2] = aux_recursive(obj.rhs,time,pred);
            mod = mod1 | mod2;
        else
            [res.rhs,mod2] = aux_recursive(obj.rhs,time,pred);
            mod = mod1 | mod2;
        end

    end  
end

% ------------------------------ END OF CODE ------------------------------
