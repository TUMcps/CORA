function [res,cert,scaling] = contains_(I,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if an STL interval contains a set or a point
%
% Syntax:
%    [res,cert,scaling] = contains_(I,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
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
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in I, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in I).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(I - I.center) + I.center contains S.
%
% Example: 
%    I1 = stlInterval(1,5,false,false);
%    I2 = stlInterval(2,4);
%
%    contains(I1,I2)
%    contains(I1,5)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% The code is not yet ready to deal with scaling or cert
cert = NaN; 
scaling = Inf;
if scalingToggle || certToggle
    throw(CORAerror('CORA:notSupported',...
        "The computation of the scaling factor or cert " + ...
        "for polynomial zonotopes is not yet implemented."));
end

% set in empty set
if representsa_(I,'emptySet',0)
    res = representsa_(S,'emptySet',0);
    return
end

if isnumeric(S)
    % point in interval
    res = (I.lower < S && S < I.upper) || ...
        (I.leftClosed && withinTol(I.lower,S,tol)) || ...
        (I.rightClosed && withinTol(I.upper,S,tol));
elseif isa(S,'stlInterval')
    % interval in interval
    upperIncluded = I.upper > S.upper || ...
        (withinTol(I.upper,S.upper,tol) && (I.rightClosed || ~S.rightClosed));
    lowerIncluded = I.lower < S.lower || ...
        (withinTol(I.lower,S.lower,tol) && (I.leftClosed || ~S.leftClosed));
    res = lowerIncluded && upperIncluded;
elseif I.leftClosed && I.rightClosed
    % we can only handle other sets exactly if I is a closed interval
    res = contains(interval(I),S,method,tol);
elseif strcmp(method,'approx')
    % for approximate containment with other sets, we simulate openness
    % by adding a small tolerance to the interval
    if ~I.leftClosed
        lower = I.lower + tol + eps;
    else
        lower = I.lower;
    end
    if ~I.rightClosed
        upper = I.upper - tol - eps;
    else
        upper = I.upper;
    end
    res = contains(interval(lower,upper),S,method,tol);
else
    throw(CORAerror('CORA:noExactAlg',I,S));
end
    

% ------------------------------ END OF CODE ------------------------------
