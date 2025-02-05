function [res,cert,scaling] = contains_(I,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if an interval contains a set or a point
%
% Syntax:
%    [res,cert,scaling] = contains_(I,S,method,maxEval,certToggle,scalingToggle,varargin)
%
% Inputs:
%    I - interval object
%    S - contSet object or single point
%    method - method used for the containment check.
%       Currently, the only available options are 'exact' and 'approx'.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of I
%       will be detected as lying in I, which can be useful to counteract
%       errors originating from floating point errors.
%    maxEval - Currently has no effect
%    certToggle - if set to 'true', cert will be computed (see below),
%       otherwise cert will be set to NaN.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below), otherwise scaling will be set to inf.
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, I is
%           guaranteed to not be contained in I, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in I).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(I - I.center) + I.center contains S.
%
% Example: 
%    I1 = interval([-1;-2],[2;3]);
%    I2 = interval([0;0],[1;1]);
%
%    contains(I1,I2)
%    contains(I1,[1;2])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, zonotope/contains_

% Authors:       Niklas Kochdumper, Mark Wetzlinger, Adrian Kulmburg
% Written:       16-May-2018
% Last update:   02-September-2019
%                15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
%                20-January-2025 (AK, added cert and scaling options)
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% init result
res = false;

% set in empty set
if representsa_(I,'emptySet',0)
    res = representsa_(S,'emptySet',0);
    cert = true;
    scaling = 0;
    return
end

% point in interval containment
if isnumeric(S)
    if scalingToggle
        c = center(I);
        scaling = abs(S)./(abs([I.inf-c I.sup-c]));
        res = scaling <= 1 + tol;
        cert = true(size(res));
    else
        res = all( (I.inf < S + tol | withinTol(I.inf,S,tol)) ...
                & (I.sup > S - tol | withinTol(I.sup,S,tol)), 1:numel(dim(I)));
        res = reshape(res,1,[]);
        cert = true(size(res));
        scaling = 0;
    end
    return

% interval in interval containment
elseif isa(S,'interval')
    % we know I is not an empty set, so only check S
    if representsa_(S,'emptySet',0)
        scaling = 0;
        cert = true;
        res = true;
        return
    end

    % compute scaling?
    if scalingToggle
        c = center(I);
        if any(isnan(c))
            nan_instances = find(isnan(c));
            for i=1:length(nan_instances)
                if I.inf(i) == -inf && I.sup(i) == inf
                    c(i) = 0;
                else
                    throw(CORAError('CORA:notDefined', ...
       "Can not determine center of the interval, as it" + ...
       "has a dimension unbounded only in one direction"));
                end
            end
        end

        I = I - c; S = S - c;
        scaling_inf = abs(S.inf./I.inf);
        scaling_sup = abs(S.sup./I.sup);

        % Need to remove NaNs, which can happen only if both
        % coordinates are the same, and are 0 or inf
        for i=1:dim(I)
            if isnan(scaling_inf)
                a=S.inf(i); b=I.inf(i);
                if b == 0
                    scaling_inf(i) = 0;
                elseif a==b
                    scaling_inf(i) = inf;
                else
                    scaling_inf(i) = 0; % This can technically never happen
                    % since I and S are both non-empty
                end
            end
            % Do the same for scaling_sup
            if isnan(scaling_sup)
                a=S.sup(i); b=I.sup(i);
                if b == 0
                    scaling_sup(i) = 0;
                elseif a==b
                    scaling_sup(i) = inf;
                else
                    scaling_sup(i) = 0;
                end
            end
        end

        scaling = max([scaling_inf; scaling_sup]);
        cert = true;
        res = scaling <= 1+tol;
        return
        
    else % do not compute scaling
        
        if all((I.sup >= S.sup | withinTol(I.sup,S.sup,tol)),"all") ...
                && all((I.inf <= S.inf | withinTol(I.inf,S.inf,tol)),"all")
            res = true;
            cert = true;
            scaling = NaN;
        else
            res = false;
            cert = true;
            scaling = NaN;
        end
        return
    end

% other set in interval containment
else

    [res,cert,scaling] = contains_(polytope(I),S,method,tol,maxEval,certToggle,scalingToggle);

end

% ------------------------------ END OF CODE ------------------------------
