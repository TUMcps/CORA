function res = modelCheckTrace(obj,x,t)
% modelCheckTrace - check if a trace satisfies an STL formula 
%
% Syntax:
%    res = modelCheckTrace(obj,x,t)
%
% Inputs:
%    obj - logic formula (class stl)
%    x - states of the trace (dimensions: [m,n])
%    t - times of the trace (dimensions: [m,1])
%
% Outputs:
%    res - formula satisfied (true) or not (false)
%
% Example: 
%    x = stl('x',2);
%    eq = aux_until(x(2) < -0.5,x(1) > 0.5,interval(0,1));
%    
%    phi = -pi/2:0.01:0;
%    x = [cos(phi'),sin(phi')];
%    t = linspace(0,1,length(phi))';
%
%    res = modelCheckTrace(eq,x,t)
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
    if abs(min(diff(t)) - max(diff(t))) > eps
        throw(CORAerror('CORA:notSupported',...
                      'Only uniformly sampled traces are supported!'));
    end

    % model check the trace using a recursive function
    vars = getVariables(obj);
    res = aux_recursive(obj,x,t,vars);
    res = res(1);
end


% Auxiliary functions -----------------------------------------------------

function res = aux_recursive(obj,x,t,vars)
% recursive function for checking traces

    if strcmp(obj.type,'&')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        rhs = aux_recursive(obj.rhs,x,t,vars);

        res = lhs & rhs;

    elseif strcmp(obj.type,'|')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        rhs = aux_recursive(obj.rhs,x,t,vars);

        res = lhs | rhs;

    elseif strcmp(obj.type,'~')

        inner = aux_recursive(obj.lhs,x,t,vars);
        res = ~inner;

    elseif strcmp(obj.type,'<')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        res = lhs < obj.rhs; 

    elseif strcmp(obj.type,'>')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        res = lhs > obj.rhs; 

    elseif strcmp(obj.type,'<=')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        res = lhs <= obj.rhs; 

    elseif strcmp(obj.type,'>=')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        res = lhs >= obj.rhs; 

    elseif strcmp(obj.type,'+')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        rhs = aux_recursive(obj.rhs,x,t,vars);

        res = lhs + rhs;

    elseif strcmp(obj.type,'-')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        rhs = aux_recursive(obj.rhs,x,t,vars);

        res = lhs - rhs;

    elseif strcmp(obj.type,'*')

        lhs = obj.lhs;
        rhs = aux_recursive(obj.rhs,x,t,vars);

        res = lhs * rhs;

    elseif strcmp(obj.type,'true')

        res = true*ones(size(t));

    elseif strcmp(obj.type,'false')

        res = false*zeros(size(t));

    elseif strcmp(obj.type,'variable')

        ind = find(contains(vars,obj.var));
        res = x(:,ind);

    elseif strcmp(obj.type,'next')

        inner = aux_recursive(obj.lhs,x,t,vars);
        ind = find(t >= obj.from);

        if ~isempty(ind)
            res = inner(ind:end);
            res = [res; false*zeros(length(t)-length(res),1)];
        else
            res = false*zeros(size(t));
        end

    elseif strcmp(obj.type,'finally')

        rhs = aux_recursive(obj.lhs,x,t,vars);
        lhs = true*ones(size(rhs));

        res = aux_until(lhs,rhs,t,obj.from,obj.to);

    elseif strcmp(obj.type,'globally')

        rhs = aux_recursive(obj.lhs,x,t,vars);
        lhs = true*ones(size(rhs));

        res = ~aux_until(lhs,~rhs,t,obj.from,obj.to);

        ind = find(t >= obj.to,1);

        if isempty(ind)
            res = false*res;
        else
            res(length(res)-ind+1:end) = false;
        end

    elseif strcmp(obj.type,'release')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        rhs = aux_recursive(obj.rhs,x,t,vars);

        res = ~aux_until(~lhs,~rhs,t,obj.from,obj.to);

        ind = find(t >= obj.to,1);

        if isempty(ind)
            res = false*res;
        else
            res(length(res)-ind+1:end) = false;
        end

    elseif strcmp(obj.type,'until')

        lhs = aux_recursive(obj.lhs,x,t,vars);
        rhs = aux_recursive(obj.rhs,x,t,vars);

        res = aux_until(lhs,rhs,t,obj.from,obj.to);

    else

        inner = aux_recursive(obj.lhs,x,t,vars);
        eval(['res = ',obj.type,'(inner);']);

    end
end

function res = aux_until(lhs,rhs,t,from,to)
% process until-operator

    ind = find(t >= from & t <= to);
    res = false*zeros(length(t),1); cnt = 1;

    while ~isempty(ind) && ind(1) <= length(t)
        
        ind = ind(ind <= length(t));

        ind1 = find(lhs(cnt:ind(end)) == false,1);
        ind2 = find(rhs(ind) == true,1);

        if ~isempty(ind2) && (isempty(ind1) || ...
                        (ind2 + ind(1) < ind1 + cnt && ind1 + cnt >= ind(1)))
            res(cnt) = true;
        end

        cnt = cnt + 1; ind = ind + 1;
    end
end

% ------------------------------ END OF CODE ------------------------------
