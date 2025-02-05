function [res,cert,scaling] = contains_(pgon1, pgon2, method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - check if polygon object pgon1 contains polygon pgon2 or a vector 
%    of points (one logical value for each point)
%
% Syntax:
%    [res,cert,scaling] = contains_(pgon1, pgon2, method,tol,maxEval,certToggle,scalingToggle,varargin)
%
% Inputs:
%    pgon1 - polygon
%    pgon2 - polygon
%    method - 'exact'
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of
%       pgon1 will be detected as lying in pgon1, which can be useful to
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
%           verified. For example, if res=false and cert=true, pgon2 is
%           guaranteed to not be contained in pgon1, whereas if res=false
%           and cert=false, nothing can be deduced (pgon2 could still be
%           contained in pgon1).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(pgon1 - pgon1.center) + pgon1.center contains pgon2.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, integration in contSet)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% The code is not yet ready to deal with scaling or cert
cert = NaN;
scaling = Inf;
if scalingToggle || certToggle
    throw(CORAerror('CORA:notSupported',...
        "The computation of the scaling factor or cert " + ...
        "for constrained polynomial zonotopes is not yet implemented."));
end

% enlarge pgon1 by tolerance
pgon1 = expandBoundaries(pgon1, tol);

if isnumeric(pgon2)

    res = isinterior(pgon1.set, pgon2(1, :), pgon2(2, :))';

elseif isa(pgon2, 'contSet')
    % convert to polygon
    pgon2 = polygon(pgon2);

    % compute union
    u = union(pgon1.set, pgon2.set);

    % check if area of pgon1 is identical to area of union
    A1 = area(pgon1.set);
    A2 = area(u);

    res = withinTol(A1, A2, tol);
else
    throw(CORAerror('CORA:notSupported', ...
        'This set representation is not supported!'));
end

end

% ------------------------------ END OF CODE ------------------------------
