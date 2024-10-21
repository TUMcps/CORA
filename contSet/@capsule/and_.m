function res = and_(C,S,varargin)
% and_ - overloads '&' operator, computes the intersection of a capsule
%    and another set or numerical vector
%
% Syntax:
%    res = C & S
%    res = and_(C,S)
%
% Inputs:
%    C - capsule object
%    S - contSet object
%
% Outputs:
%    res - intersection
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Mark Wetzlinger
% Written:       28-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call function with lower precedence
if isa(S,'contSet') && S.precedence < C.precedence
    res = and_(S,C,varargin{:});
    return
end

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",C,S));

% ------------------------------ END OF CODE ------------------------------
