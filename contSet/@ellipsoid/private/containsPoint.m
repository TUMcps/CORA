function [B,Val] = containsPoint(E,Y)
% containsPoint - gives an array of boolean values indiciating whether
%    points Y are contained in the ellipsoid
%
% Syntax:  
%    [B,Val] = containsPoint(E,Y) gives an array of boolean values indiciating
%     whether points Y are contained in the ellipsoid
%
% Inputs:
%    E - ellipsoids object
%    Y - Points
%
% Outputs:
%    B - boolean values indiciating whether
%        points Y are contained in the ellipsoid
%    Val-if contained, Value indicates the relative distance to the center of E: Val<=1
%    <=> contained (=1: on boundary); otherwise: inf
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
% Last update:  16-March-2021 (complete rewrite)
% Last revision:---

%------------- BEGIN CODE --------------
n = dim(E);
if ~isa(Y,'double')
    error('Second argument must be a double matrix');
end
[m,N] = size(Y);
if m~=n
    error('First dimension of second input does not match ellipsoid dimension!');
end


B = false(1,N);
Val = inf(1,N);
ind_rem_eq = true(1,N);
if E.isdegenerate
    [T,~,~] = svd(E.Q);
    E = T'*E;
    Y = T'*Y;
    % save remainder
    x_rem = E.q(E.rank+1:end);
    Y_rem = Y(E.rank+1:end,:);
    % check whether x_rem==Y_rem (those that do not fullfill that are
    % already not contained)
    % indices of B which might be contained
    ind_rem_eq = withinTol(Y_rem,repmat(x_rem,1,size(Y_rem,2)),E.TOL);
    % if only center remains
    if rank(E)==0
        B(ind_rem_eq) = true;
        Val(ind_rem_eq) = 1;
        return;
    end
    % project so that E is no longer degenerate
    E = project(E,1:E.rank);
    Y = Y(1:E.rank,:);
end
% convert mask to indices
tmp = 1:N;
ii_eq_rem = tmp(ind_rem_eq);
% now, E is fulldimensional
for i=ii_eq_rem
    % simply check using ellipsoid equation
    val_i = (Y(:,i)-E.q)'*inv(E.Q)*(Y(:,i)-E.q);
    B(i) = val_i <= 1+E.TOL;
    Val(i) = B(i)*val_i;
end
%------------- END OF CODE --------------