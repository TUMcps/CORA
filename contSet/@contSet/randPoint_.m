function p = randPoint_(varargin)
% randPoint_ - generates a random point within a given continuous set
%    (internal use, see also contSet/randPoint)
%
% Syntax:
%    p = randPoint_(S, N)
%    p = randPoint_(S, N, type)
%    p = randPoint_(S, 'all', 'extreme')
%
% Inputs:
%    pZ - polyZonotope object
%    N - number of random points (default: 1)
%    type - type of the random point ('standard' or 'extreme')
%           Note that for 'extreme', the generated points might not all be
%           extremal for all sets
%
% Outputs:
%    p - random points in R^n
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
