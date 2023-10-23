function res = isYalmipInstalled()
% isYalmipInstalled - checks if YALMIP [1] is installed
%
% Syntax:
%    res = isYalmipInstalled
%
% Inputs:
%    -
%
% Outputs:
%    res - true if installed, false otherwise
%
% Example: 
%    -
%
% References:
%    [1] Löfberg, J., 2004, September. YALMIP: A toolbox for 
%        modeling and optimization in MATLAB. In Proceedings of 
%        the CACSD Conference (Vol. 3).
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @zonotope/minnorm.m

% Authors:       Victor Gassmann
% Written:       15-October-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    sdpvar;
    res = true;
catch
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
