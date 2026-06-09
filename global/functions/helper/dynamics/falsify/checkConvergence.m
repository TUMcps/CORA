function res = checkConvergence(y)
% checkConvergence - check if a series of values converges to a value 
%    greater 0 (res = 1) or not (res = 0)
%
% Syntax:
%    res = checkConvergence(y)
%
% Inputs:
%    y - vector representing a series of values
%
% Outputs:
%    res - true if converging, false if not
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       20-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    res = 0;

    % check if there are enough values to judge convergence
    if length(y) < 3
        return;
    end

    % check if the value is increasing
    d = diff(y);

    if any(d > 0)
        res = 1; return;
    end

    % check if last point is an outlier
    d = diff(y);
    a = d(2:end)./d(1:end-1);

    if a(end) > 1
        return;
    end

    % estimate limit values
    l = zeros(length(y)-2,1);

    for i = 3:length(y)
        l(i-2) = aux_estimateLimitValue(y(1:i));
    end

    d = diff(l);

    if any(d >= 0) || (l(end) > 0 && any(l < 0))
        res = 1;
    end
end


% Auxiliary functions -----------------------------------------------------

function l = aux_estimateLimitValue(y)
% estimate the limit value of a sequence of decreasing values y based on 
% the geometric series

    % too less values to do a meaningful estimation
    if length(y) < 3
        l = 0; return;
    end

    % catch the case where the values are increasing
    d = diff(y);

    if any(d > 0)
        l = min(y); return;
    end

    % estimate limit value based on the "geometric series"
    a = d(2:end)./d(1:end-1);
    
    if all(a > 1-eps)
        l = 0; return;
    end

    a = a(a < 1-eps);
    s = (1:length(a))';
    a = a.*(s/sum(s));
    a = sum(a);

    l = y(end-1) + d(end)/(1-a); 
end

% ------------------------------ END OF CODE ------------------------------
