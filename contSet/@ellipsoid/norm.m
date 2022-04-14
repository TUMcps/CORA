function val = norm(E,type)
% norm - compute the maximum Euclidean norm of an ellipsoid
%
% Syntax:  
%    val = norm(E,type)
%
% Inputs:
%    E    - Ellipsoid object 
%    type - norm type (2)
%
% Outputs:
%    val - value of the maximum norm
%
% Example: 
%    E = ellipsoid.generateRandom(false,2);
%    val = norm(ellipsoid(E.Q,zeros(size(E.q))));
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      20-November-2019
% Last update:  31-July-2020
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('type','var')
    type = 2;
end
if type~=2
    error('Only implemented for Euclidean norm');
end
if ~all(E.q==0)
    error('Not implemented for non-zero center yet');
end
% transform into eigenspace
lmax = max(eig(E.Q));
val = sqrt(lmax);
%------------- END OF CODE --------------