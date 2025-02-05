function [res,cert,scaling] = contains_(O,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if an empty set contains a set or a point
%
% Syntax:
%    [res,cert,scaling] = contains_(O,S,method,tol,maxEval,certToggle,scalingToggle)
%
% Inputs:
%    O - emptySet object
%    S - contSet object or numerical vector
%    method - method used for the containment check.
%       Currently, the only available options are 'exact' and 'approx'.
%    tol - tolerance for the containment check; not used here.
%    maxEval - Currently has no effect
%    certToggle - if set to 'true', cert will be computed (see below),
%       otherwise cert will be set to NaN.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below), otherwise scaling will be set to inf.
%
% Outputs:
%    res - true/false
%
% Example: 
%    O = emptySet(2);
%    p = [1;1];
%    [res,cert,scaling] = contains(O,p);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (MW, rename contains_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dimensions are already checked...
res = false;
cert = true;
scaling = inf;
if isa(S,'emptySet')
    % empty set contains the empty set
    res = true;
    scaling = 0;

elseif isa(S,'contSet') && representsa(S,'emptySet',eps)
    % empty set contains contSet objects if they also represent the empty
    % set
    res = true;
    scaling = 0;

elseif isnumeric(S) && isempty(S)
    % empty set contains empty vectors
    res = true;
    scaling = 0;
end

% ------------------------------ END OF CODE ------------------------------
