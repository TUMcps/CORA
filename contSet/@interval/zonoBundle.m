function zB = zonoBundle(I)
% zonoBundle - converts an interval into a zonotope bundle
%
% Syntax:
%    zB = zonoBundle(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    zB - zonoBundle object
%
% Example:
%    I = interval([-2;-3;-4],[4;3;2]);
%    zB = zonoBundle(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/zonoBundle

% Authors:       Niklas Kochdumper
% Written:       26-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

zB = zonoBundle({zonotope(I)});

% ------------------------------ END OF CODE ------------------------------
