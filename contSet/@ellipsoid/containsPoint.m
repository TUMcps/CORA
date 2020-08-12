function [B] = containsPoint(E,Y)
% containsPoint - gives an array of boolean values indiciating whether
%    points Y are contained in the ellipsoid
%
% Syntax:  
%    B = containsPoint(E,Y) gives an array of boolean values indiciating
%     whether points Y are contained in the ellipsoid
%
% Inputs:
%    E - ellipsoids object
%    Y - Points
%
% Outputs:
%    B - boolean values indiciating whether
%        points Y are contained in the ellipsoid
%
% Example: 
%    t = linspace(0,2*pi,1000);
%    Y = [cos(t);sin(t)];
%    E = ellipsoid([1,0;0,1/2],[1;1]);
%    B = containsPoint(E,Y);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%If dimension is lower, we use the reduced model
B = false(size(Y,2),1);
[V,D] = eig(E.Q);
[~,ind] = sort(diag(D),'descend');
D = D(ind,ind);
V = V(:,ind);

if E.dim<length(E.Q)
    V_n = V(:,E.dim+1:end);
    V = V(:,1:E.dim);
    D=D(1:E.dim,1:E.dim);
end

for i=1:length(B)
    y = Y(:,i);
    y0 = y-E.q;
    if E.dim<length(E.Q) && ~all(abs(V_n'*y0)<=E.TOL)
        continue;
    end
    y_t = V'*y0;
    B(i) = y_t'*diag(1./diag(D))*y_t<=1+E.TOL;
end

%------------- END OF CODE --------------