function r = radius(obj, varargin)
% radius - Computes radius of an enclosing hyperball of an ellipsoid
%
% Syntax:  
%    r = radius(obj) returns the largest radius
%    r = radius(obj,i) returns the i largest radii
%
% Inputs:
%    obj - ellipsoid object
%
% Outputs:
%    r - radius of enclosing hyperball/vector of radii
%
% Example: 
%    E = ellipsoid([1,0.5;0.5,3],[1;-1]);
%    r = radius(E);
%
%    figure
%    hold on
%    plot(E,[1,2],'r');
%    Ecirc = ellipsoid(r^2 * eye(2),center(E));
%    plot(Ecirc,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      05-March-2021
% Last update:  19-March-2021 (VG: empty case added)
% Last revision:---

%------------- BEGIN CODE --------------

% no further input
if nargin==1
    i = 1; % only largest radius considered
elseif nargin==2
    i = varargin{1};
else
    disp('too many inputs');
end

if isempty(obj)
    r = [];
    return;
end

% compute eigenvalues
d = eigs(obj.Q,i); % since we use Q^{-1} as a shape matrix

%compute radius
r = sqrt(d); % since we use Q^{-1} as a shape matrix

%------------- END OF CODE --------------