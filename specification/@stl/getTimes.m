function [timeInt,timePoint] = getTimes(obj)
% getTimes - get times at which the predicates appear in the formula
%
% Syntax:
%    [timeInt,timePoint] = getTimes(obj)
%
% Inputs:
%    obj - logic formula (class stl)
%
% Outputs:
%    timeInt - list storing for each predicates the time intervals in which
%              the predicates appears in the formula
%    timePoint - list storing for each predicates the time points in which
%                the predicates appears in the formula
%
% Example: 
%    x = stl('x',2);
%    eq = x(1) < 5 | globally(x(2) < 3,interval(1,2));
%    eq = assignIdentifiers(eq);
%    [timeInt,timePoint] = getTimes(eq);
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
    
    % determine times with a recursive function
    ind = aux_recursive(obj,0);

    timePoint = {}; timeInt = {};

    if isfield(ind,'timePoint')
        timePoint = ind.timePoint;
    end

    if isfield(ind,'timeInt')
        timeInt = ind.timeInt;
    end

    % bring lists to a common format
    if length(timeInt) > length(timePoint)
        timePoint = [timePoint; cell(length(timeInt)-length(timePoint),1)];
    elseif length(timePoint) > length(timeInt)
        timeInt = [timeInt; cell(length(timePoint)-length(timeInt),1)];
    end
end


% Auxiliary functions -----------------------------------------------------

function ind = aux_recursive(obj,time)
% recursive function to determine times where the predicates are active

    if isempty(obj.id) && (strcmp(obj.type,'&') || strcmp(obj.type,'|'))

        ind1 = aux_recursive(obj.lhs,time);
        ind2 = aux_recursive(obj.rhs,time);
        
        ind = aux_mergeIndices(ind1,ind2);

    elseif ~obj.temporal

        if ~isempty(obj.id)
            if isnumeric(time)
                ind.timePoint{obj.id,1} = {time};
            else
                ind.timeInt{obj.id,1} = {time};
            end
        else
            ind = [];
        end

    elseif strcmp(obj.type,'next')
        
        ind = aux_recursive(obj.lhs,obj.from+time);

    elseif strcmp(obj.type,'globally') || strcmp(obj.type,'finally')

        ind = aux_recursive(obj.lhs,interval(obj.from,obj.to) + time);

    elseif strcmp(obj.type,'until') || strcmp(obj.type,'release') 

        ind1 = aux_recursive(obj.lhs,interval(obj.from,obj.to) + time);
        ind2 = aux_recursive(obj.rhs,interval(obj.from,obj.to) + time);
        
        ind = aux_mergeIndices(ind1,ind2);
    end  
end

function ind = aux_mergeIndices(ind1,ind2)
% merge the lists of times

    if isempty(ind1)
        ind = ind2; return;
    elseif isempty(ind2)
        ind = ind1; return;
    end

    if ~isfield(ind1,'timePoint') && isfield(ind2,'timePoint')
        ind.timePoint = ind2.timePoint;
    elseif ~isfield(ind2,'timePoint') && isfield(ind1,'timePoint')
        ind.timePoint = ind1.timePoint;
    elseif isfield(ind1,'timePoint') && isfield(ind2,'timePoint')
        l1 = length(ind1.timePoint); l2 = length(ind2.timePoint);
        ind.timePoint = cell(max(l1,l2),1);
        for i = 1:l1
            ind.timePoint{i} = ind1.timePoint{i};
        end
        for i = 1:l2
            ind.timePoint{i} = [ind.timePoint{i},ind2.timePoint{i}];
        end
    end

    if ~isfield(ind1,'timeInt') && isfield(ind2,'timeInt')
        ind.timeInt = ind2.timeInt;
    elseif ~isfield(ind2,'timeInt') && isfield(ind1,'timeInt')
        ind.timeInt = ind1.timeInt;
    elseif isfield(ind1,'timeInt') && isfield(ind2,'timeInt')
        l1 = length(ind1.timeInt); l2 = length(ind2.timeInt);
        ind.timeInt = cell(max(l1,l2),1);
        for i = 1:l1
            ind.timeInt{i} = ind1.timeInt{i};
        end
        for i = 1:l2
            ind.timeInt{i} = [ind.timeInt{i},ind2.timeInt{i}];
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
