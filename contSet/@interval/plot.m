function han = plot(I,varargin)
% plot - plots a projection of an interval 
%
% Syntax:
%    han = plot(I)
%    han = plot(I,dims)
%    han = plot(I,dims,type)
%
% Inputs:
%    I - interval object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
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

%plot interval via conversion to zonotope
if nargin == 1
    han = plot(zonotope(I));
else
    dims = varargin{1};
    % we need to project the interval first since intervals with -Inf/Inf
    % values exist which result in zonotopes with NaN centers...
    Z = zonotope(project(I,dims));

    if length(dims) == 1
        han = plot(Z,1,varargin{2:end});
    elseif length(dims) == 2
        han = plot(Z,[1,2],varargin{2:end});
    elseif length(dims) == 3
        han = plot(Z,[1,2,3],varargin{2:end});
    end
end

if nargout == 0
    clear han;
end

%------------- END OF CODE --------------