function [Y,L_t] = boundary(E,N)
% boundary - Computes the boundary points for all directions contained in L
%
% Syntax:  
%    S = Y(E,N) Computes N samples of boundary points of E
%
% Inputs:
%    E - ellipsoid object
%    N - #Samples
%
% Outputs:
%    Y - boundary points
%    L - Sampled directions (unit vectors)
%
% Example: 
%    t = linspace(0,2*pi,1000);
%    E = ellipsoid([1,0;0,1/2],[1;1]);
%    Y = boundary(E,1000);
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
%NOTICE: Does NOT generate uniformly distributed boundary points
if N==0
    error('Please provide N>0');
end

if length(E.Q)==1
    L_t = [];
    %If we only have a scalar unknown (i.e. length(Obj.Q)==1): If
    %Obj.Q==0, then a point x in the ellipsoid has to satisfy l*x<=l*q
    %for every l in R. Therefore x=q.
    if E.Q>E.TOL
        Y = linspace(-sqrt(E.Q)+E.q,sqrt(E.Q)+E.q,N);
    else
        disp('Returning N times the only boundary point for Q=0');
        Y = repmat(E.q,1,N);
    end
    return;
end

[V,D] = eig(E.Q);
[~,ind] = sort(diag(D),'descend');
D = D(ind,ind);
V = V(:,ind);
if E.dim<length(E.Q)
    V = V(:,1:E.dim);
    D = D(1:E.dim,1:E.dim);
end
X_t = randn(length(D),N);
%uniformly sample the n-hypersphere
L_t = 1./sqrt(sum(X_t.^2,1)).*X_t;
Y = zeros(length(E.Q));

for i=1:size(L_t,2)
    l_t = L_t(:,i);
    y_t = D*l_t/sqrt(l_t'*D*l_t);
    Y(:,i) = V*y_t+E.q;
end

%------------- END OF CODE --------------