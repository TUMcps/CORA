function res = sign(I)
% sign - overloaded built-in signum function; an exact evaluation would
%    only contain the values -1, 0, and 1, thus, our result represents an
%    outer approximation
%
% Syntax:
%    res = sign(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example:
%    I = interval([-2;3],[3;4]);
%    res = sign(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Sebastian Mair, Mark Wetzlinger
% Written:       17-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = interval(sign(I.inf), sign(I.sup));

% ------------------------------ END OF CODE ------------------------------
