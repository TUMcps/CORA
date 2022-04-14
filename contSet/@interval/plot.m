function han = plot(obj,varargin)
% plot - Plots 2-dimensional projection of an interval 
%
% Syntax:
%    plot(obj)
%    plot(obj,dims)
%    plot(obj,dims,'r',...)
%
% Inputs:
%    obj - interval object
%    dims - (optional) dimensions that should be projected
%    type - (optional) plot type (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the plotted object
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
Z = zonotope(obj);

%plot zonotope
if nargin == 1
    han = plot(Z);
else
    han = plot(Z,varargin{:});
end

%------------- END OF CODE --------------