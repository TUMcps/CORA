function E = orEllipsoidOA(E,E_cell)
% minus - Computes the outer approximation of the union between
%           ellipsoids
%
% Syntax:  
%    E = unionEllipsoidOA(E,E_cell)
%
% Inputs:
%    E1 - ellipsoid object
%    E_cell - ellipsoid array
%
% Outputs:
%    E - ellipsoid after union
%
% References:
%   [1] Boyd et al. Convex Optimization
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      15-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% collapse cell array
E_cell = [{E};E_cell];
N = length(E_cell);
n = dim(E);

% normalize to prevent numerical issues
max_val = max(cellfun(@(cc)max(svd(cc.Q)),E_cell));


fac = 0.001;
th = fac*max_val;
th(th==0) = fac;
for i=1:N
    Ei = E_cell{i};
    if Ei.isdegenerate
        % if Ei is degenerate, add small pertubation (possible since we
        % compute an overapproximation)
        nd_i = rank(Ei);
        [Ti,Si,~] = svd(Ei.Q);
        si = diag(Si);
        % choose singular value such that rcond not too small
        % bloat to remove degeneracy
        % use max values from other ellipsoids to avoid numerical problems
        % when solving the problem below
        Si = diag([si(1:nd_i);th*ones(n-nd_i,1)]);
        E_cell{i} = ellipsoid(Ti*Si*Ti',Ei.q);
    end
end

% find minimum volume ellipsoid spanning union [1]
A2 = sdpvar(n); % pd-ness of A ensured by yalmip internally
bt = sdpvar(n,1);
l = sdpvar(N,1);
f_obj = -1/2*geomean(A2);
C = [];
for i=1:N
    Qiinv = inv(E_cell{i}.Q);
    qi = E_cell{i}.q;
    Ci = [A2-l(i)*Qiinv,    bt+l(i)*Qiinv*qi,            zeros(n);
          (bt+l(i)*Qiinv*qi)',-1-l(i)*(qi'*Qiinv*qi-1),  bt';
          zeros(n),          bt,                         -A2];
    C = [C, Ci<=0];
end
opts = sdpsettings;
opts.verbose = 0;
if exist('sdpt3','file')
    opts.solver = 'sdpt3';
end
sol = optimize([C,l>=0],f_obj,opts);
if any(sol.problem == [-2,-3,-4])
    error('This function requires an SDP solver compatible with YALMIP!');
elseif sol.problem ~= 0
    throw(errOptNumIssue());
end

% extract ellipsoid
Q = inv(value(A2));
A = sqrtm(value(A2));
b = A\value(bt);
q = -Q*A*b;

% backtransform
E = ellipsoid(Q,q);
%------------- END OF CODE --------------