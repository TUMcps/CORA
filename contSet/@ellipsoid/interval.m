function I = interval(E)
% interval - Overapproximates an ellipsoid by an interval hull
%
% Syntax:  
%    I = interval(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    I - interval object
%
% Example: 
%    E = ellipsoid.generateRandom(0);
%    I = interval(E);
%    plot(E);
%    hold on
%    plot(I);
%
% Other m-files required: interval (zonotope)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
n = E.dim;
E0 = ellipsoid(E.Q,zeros(size(E.q)));
dI = zeros(n,1);
Idty = eye(n);
% compute the width of the ellipsoid in each dimension 
% using the support function
for i=1:n
    dI(i) = supportFunc(E0,Idty(:,i));
end
% construct the resulting interval
I = interval(-dI,dI) + E.q;
%------------- END OF CODE --------------