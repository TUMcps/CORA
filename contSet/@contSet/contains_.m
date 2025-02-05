function [res,cert,scaling] = contains_(S1,S2,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if a set contains another set or a point
%    (internal use, see also contSet/contains)
%
% Syntax:
%    [res,cert,scaling] = contains_(S1,S2,method,tol,maxEval,certToggle,scalingToggle)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object or numeric
%    method - method used for the containment check.
%       Currently, the only available options are 'exact' and 'approx'.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of
%       S1 will be detected as lying in S1, which can be useful to
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
%           verified. For example, if res=false and cert=true, S2 is
%           guaranteed to not be contained in S1, whereas if res=false and
%           cert=false, nothing can be deduced (S2 could still be
%           contained in S1).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(S1 - S1.center) + S1.center contains S2.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
