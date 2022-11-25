function plot(varargin)
% plot - Plots 2-dimensional projection of an interval matrix
%
% Syntax:  
%    plot(obj,dimensions)
%
% Inputs:
%    obj - interval matrix
%    dimensions - dimensions that should be projected (optional) 
%
% Outputs:
%    none
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      22-June-2010
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%convert to interval matrix to interval
IH=interval(varargin{1});
    
%plot interval
if nargin==1
    plot(IH);
elseif nargin>=2
    plot(IH,varargin{2:end});
end

%------------- END OF CODE --------------