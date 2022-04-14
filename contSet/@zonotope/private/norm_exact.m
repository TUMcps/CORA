function [val,x] = norm_exact(Z,type)
% diffnorm - computes the exact maximum norm
%
% Syntax:  
%    [val,x] = norm_exact(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%   val - norm value of vertex with biggest distance from the center
%   x   - vertex attaining maximum norm
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: minnorm

% Author:       Victor Gassmann
% Written:      18-September-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~isYalmipInstalled()
    error('YALMIP must be on the MATLAB search path to use this function');
elseif str2double(yalmip('version'))<20190425 % version: 25.04.2020
    error('YALMIP version >=20190425 required');
end

if ~exist('type','var')
    type = 2;
end
if type~=2
    error('Only euclidean norm implemented so far');
end
%%%ATTENTION: (1) is a Binary Quadratic Program, thus not very
%%%scalable!
%(1) norm(Z)^2 = max_{u\in{-1,1}^(.)} (c+G*u)'*(c+G*u) 
G = generators(Z);
if isempty(G)
    x = center(Z);
    val = sqrt(x'*x);
    return;
end
[~,m] = size(G);
GG = G'*G;
c = center(Z);
%equivalent transformation of (1) to -min_{u\in{-1,1}^m} u'*M*u - 2*c'*G*u - lmax*m
lmax = max(eig(GG));
M = lmax*eye(m) - GG;
b = binvar(m,1);
obj = 4*(b-0.5)'*M*(b-0.5) - 4*c'*G*(b-0.5);
options = sdpsettings('verbose',0);
%use gurobi for MUCH faster solve times
%options.solver = 'gurobi';
options.solver = 'bnb';
options.bnb.solver = 'quadprog';
options.bnb.maxiter = Inf;
options.bnb.maxtime = 3600;
options.gaptol = 1e-10;
%solve optimization problem
optimize([], obj, options);
val = sqrt(value(-obj + m*lmax));
x = G*value(2*(b-0.5))+c;
%------------- END OF CODE --------------