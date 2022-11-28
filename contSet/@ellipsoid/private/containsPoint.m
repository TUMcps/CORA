function [B,val] = containsPoint(E,Y)
% containsPoint - gives an array of boolean values indicating whether the
%    points Y are contained in the ellipsoid
%
% Syntax:
%    [B,val] = containsPoint(E,Y) 
%
% Inputs:
%    E - ellipsoid object
%    Y - points
%
% Outputs:
%    B - true/false indicating whether points are contained in the ellipsoid
%    val - robustness value of each point with, i.e., gives the relative
%          distance to the center of E: val <= 1 <=> contained
%          (val = 1 => on boundary); otherwise: Inf
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

% read dimension
n = dim(E);

if ~isa(Y,'double')
    throw(CORAerror('CORA:wrongValue','second','be a double matrix'));
end

[m,N] = size(Y);
if m~=n
    throw(CORAerror('CORA:dimensionMismatch',E,Y));
end

% init values
B = false(1,N);
val = inf(1,N);
ind_rem_eq = true(1,N);

if ~isFullDim(E)
    [T,~,~] = svd(E.Q);
    E = T'*E;
    Y = T'*Y;
    % save remainder
    rankE = rank(E);
    x_rem = E.q(rankE+1:end);
    Y_rem = Y(rankE+1:end,:);
    % check whether x_rem==Y_rem (those that do not fullfill that are
    % already not contained)
    % indices of B which might be contained
    ind_rem_eq = withinTol(Y_rem,repmat(x_rem,1,size(Y_rem,2)),E.TOL);
    % if only center remains
    if rankE==0
        B(ind_rem_eq) = true;
        val(ind_rem_eq) = 1;
        return;
    end
    % project so that E is no longer degenerate
    E = project(E,1:rankE);
    Y = Y(1:rankE,:);
end

% convert mask to indices
tmp = 1:N;
ii_eq_rem = tmp(ind_rem_eq);

% now, E is fulldimensional
for i=ii_eq_rem
    % simply check using ellipsoid equation
    val_i = (Y(:,i)-E.q)'*inv(E.Q)*(Y(:,i)-E.q);
    B(i) = val_i <= 1+E.TOL;
    val(i) = B(i)*val_i;
end

%------------- END OF CODE --------------