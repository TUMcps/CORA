function [A,b,Ae,be] = priv_plus_minus_vector(A,b,Ae,be,v)
% priv_plus_minus_vector - computes the translation of a polytope by a
%    vector, assumed 'P+v' (for 'P-v', simply hand over -v)
%
% Syntax:
%    [A,b,Ae,be] = priv_plus_minus_vector(A,b,Ae,be,v)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    v - numeric
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute shift
b = b + A*v;
be = be + Ae*v;

% ------------------------------ END OF CODE ------------------------------
