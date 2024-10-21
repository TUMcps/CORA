function I = interval(obj,varargin)
% interval - overapproximates an STL interval as a closed interval
%
% Syntax:
%    I = interval(obj)
%
% Inputs:
%    obj - stlInterval object
%
% Outputs:
%    I - interval object
%
% Example:
%    int = stlInterval(1,2,false,false);
%    I = interval(int);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

I = interval(obj.lower,obj.upper);

% ------------------------------ END OF CODE ------------------------------
