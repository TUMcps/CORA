function P_out = normalizeConstraints(P,varargin)
% normalizeConstraints - normalizes all inequality/equality constraints;
%    type 'b' (default) normalizes the entries in the offset vector b, so
%    that we have
%       A(i,:) x <= 1
%       A(i,:) x <= 0
%       A(i,:) x <= -1
%    for the inequality constraints and
%       Ae(i,:) x = 1
%       Ae(i,:) x = 0
%    for the equality constraints; alternatively, using type 'A', the norm
%    of the constraints in A and Ae is normalized to 1.
%
% Syntax:
%    P = normalizeConstraints(P)
%    P = normalizeConstraints(P,type)
%
% Inputs:
%    P - polytope object
%    type - (optional) 'b', 'be' (default): normalize offset vectors b and be
%                      'A', 'Ae': normalize norm of constraints in A and Ae to 1
%
% Outputs:
%    P_out - polytope object
%
% Example: 
%    A = [2 1; -2 3; -2 0; -1 -4; 2 -3; 4 -1];
%    b = [2 3 -1 2 0 1]';
%    P = polytope(A,b);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       01-December-2022
% Last update:   13-December-2022 (MW, add A normalization)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% can only normalize constraints if constraints are given
if ~P.isHRep.val
    throw(CORAerror('CORA:specialError',...
        'Halfspace representation not computed... nothing to normalize!'));
end

% set default type
type = setDefaultValues({'b'},varargin);

% check input arguments
inputArgsCheck({{P,'att','polytope'}, ...
                {type,'str',{'b','be','A','Ae'}}});

% copy polytope including properties
P_out = polytope(P);

% normalize constraints
[P_out.A_.val,P_out.b_.val,P_out.Ae_.val,P_out.be_.val] = ...
    priv_normalizeConstraints(P.A,P.b,P.Ae,P.be,type);

% ------------------------------ END OF CODE ------------------------------
