function [res,cert,scaling] = contains_(ls,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if a level set contains a set or a point
%
% Syntax:
%    [res,cert,scaling] = contains_(ls,S,method,tol,maxEval,certToggle,scalingToggle)
%
% Inputs:
%    ls - levelSet object
%    S - contSet object or single point
%    method - method used for the containment check.
%       Currently, the only available options are 'exact' and 'approx'.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of
%       ls will be detected as lying in ls, which can be useful to
%       counteract errors originating from floating point errors.
%    maxEval - Currently has no effect
%    certToggle - if set to 'true', cert will be computed (see below),
%       otherwise cert will be set to NaN.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below), otherwise scaling will be set to inf.
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in ls, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in ls).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(ls - ls.center) + ls.center contains S.
%
% Example: 
%    syms x y
%    eq = sin(x) + y;
%    ls = levelSet(eq,[x;y],'<=');
%
%    I1 = interval([0.7;-0.3],[1.3;0.3]);
%    I2 = interval([-1.3;-0.3],[-0.7;0.3]);
%
%    contains(ls,I1)
%    contains(ls,I2)
%
%    figure; hold on; xlim([-1.5,1.5]); ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(I1,[1,2],'FaceColor','r');
%
%    figure; hold on; xlim([-1.5,1.5]); ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(I2,[1,2],'FaceColor','g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, polytope/contains_

% Authors:       Niklas Kochdumper
% Written:       19-July-2019
% Last update:   15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

if representsa(S, 'emptySet', tol)
    % Empty set is always contained
    res = true;
    cert = true;
    scaling = 0;
    return
elseif representsa(S, 'fullspace', tol)
    % Fullspace is never contained, since a cPZ is compact
    res = false;
    cert = true;
    scaling = inf;
    return
end

% The code is not yet ready to deal with scaling
cert = false;
scaling = Inf;
if scalingToggle
    throw(CORAerror('CORA:notSupported',...
        "The computation of the scaling factor or cert " + ...
        "for constrained polynomial zonotopes is not yet implemented."));
end

% set or single point 
if ~isnumeric(S)                                     % set

    % check type of level set
    if strcmp(ls.compOp,'==')
        throw(CORAerror('CORA:noops',ls,S));
    end

    % interval over-approximation
    I = interval(S); 

    % evaluate non-linear function with interval arithmetic
    I = ls.funHan(I);

    ub = supremum(I);

    % multiple inequality constraints or not
    if ~iscell(ls.compOp)

        % switch case for the different types of level sets
        if strcmp(ls.compOp,'<=')
            res = ub < 0 | withinTol(ub,0,tol);
            cert = res;
        else
            res = ub < 0;
            cert = res;
        end

    else

        resVec = false(length(ls.compOp),1);

        % loop over all inequality constraints
        for i = 1:length(ls.compOp)

            if strcmp(ls.compOp{i},'<=')
                resVec(i) = ub(i) < 0 | withinTol(ub(i),0);
            else
                resVec(i) = ub(i) < 0;
            end
        end

        res = all(resVec);
        cert = res;
    end

else                                                    % array of points
    cert = true;
    % init resulting logical array
    res = false(1,size(S,2));

    for j=1:size(S,2)
    
        % evaluate nonlinear function
        val = ls.funHan(S);
        
         % multiple inequality constraints or not
        if ~iscell(ls.compOp)
    
            % switch case for the different types of level sets
            if strcmp(ls.compOp,'==')
                tmp = abs(val);
                res(j) = all(tmp < tol | withinTol(tmp,tol));
            elseif strcmp(ls.compOp,'<=')
                res(j) = all(val < tol | withinTol(val,tol));
            else
                res(j) = all(val < tol);
            end
    
        else
    
            resVec = false(length(ls.compOp),1);
    
            % loop over all inequality constraints
            for i = 1:length(ls.compOp)
    
                if strcmp(ls.compOp{i},'==')
                    tmp = abs(val(i));
                    resVec(i) = tmp < tol | withinTol(tmp,tol);
                elseif strcmp(ls.compOp{i},'<=')
                    resVec(i) = val(i) < tol | withinTol(val(i),tol);
                else
                    resVec(i) = val(i) < tol;
                end
            end
    
            res(j) = all(resVec);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
