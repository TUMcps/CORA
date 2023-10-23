function rem = getRem(tay)
% getRem - returns remainder
%
% Syntax:
%    rem = getRem(tay)
%
% Inputs:
%    tay - a Taylor model
%
% Outputs:
%    rem - an interval 
%
% Example:
%    syms x y
%    func = sin(x+y) + exp(-x) + x*y;
%    tay = taylm(func,interval([1;3],[2;4]),4);
%    getRem(tay)
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Dmitry Grebenyuk
% Written:       06-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rem = tay.remainder;

% ------------------------------ END OF CODE ------------------------------
