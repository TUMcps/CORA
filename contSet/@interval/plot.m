function plot(varargin)
% plot - Plots 2-dimensional projection of an interval 
%
% Syntax:  
%    plot(obj,dimensions)
%
% Inputs:
%    obj - interval object
%    dimensions - dimensions that should be projected (optional) 
%    type - plot type
%
% Outputs:
%    none
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    plot(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      31-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%convert to zonotope
Z = zonotope(varargin{1});
    
%plot zonotope
if nargin == 1
    plot(Z);
else
    plot(Z,varargin{2:end});
end

%------------- END OF CODE --------------