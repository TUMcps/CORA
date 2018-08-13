function [Z] = zonotope(varargin)
% mSigma - converts a probabilistic zonotope to a common zonotope where for
% each generator, a m-sigma interval is taken
%
% Syntax:  
%    [Z] = zonotope(obj,m)
%
% Inputs:
%    obj - probabilistic zonotope object
%    m - integer determining the size of the intervals
%
% Outputs:
%    Z - standard zonotope
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 03-August-2007 
% Last update: 26-February-2008
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin==1
    obj=varargin{1};
    m=obj.gamma;
elseif nargin==2
    obj=varargin{1};
    m=varargin{2};   
end

%reduce probabilistic zonotope
[pZred]=probReduce(obj);

%new generators
newG=m*pZred.g;

%build new zonotope
c=obj.Z(:,1);
G=[obj.Z(:,2:end),newG];
Z=zonotope([c,G]);

%------------- END OF CODE --------------