function P_out = reduceOverDomain(P,dom,varargin)
% reduceOverDomain - reduce the number of inequality constraints in a
%    polytope to a single inequality constraint that either outer- or
%    inner-approximates the original polytope within a given domain;
%    caution: this is not equivalent to a rigorous reduction over the R^n!
%    (note: equality constraints not supported)
%
% Syntax:
%    hs = reduceOverDomain(P,dom)
%    hs = reduceOverDomain(P,dom,type)
%
% Inputs:
%    P - polytope
%    dom - domain of values (class contSet)
%    type - inner-approximation ("inner") or outer-approximation ("outer")
%
% Outputs:
%    P_out - reduced polytope
% 
% Example: 
%    dom = interval(-ones(2,1),ones(2,1));
%    P = polytope([-1 -1; -1.2 -0.9; -0.9 -1.1],[-1;-1;-1]);
%
%    P_inner = reduceOverDomain(P,dom,'inner');
%    P_outer = reduceOverDomain(P,dom,'outer');
%
%    figure; hold on; box on;
%    xlim([-2,2]); ylim([-2,2]);
%    plot(dom);
%    plot(P,[1,2],'r','FaceAlpha',0.5);
%    plot(polytope([],[],P_inner.A,P_inner.b),[1,2],'g');
%    plot(polytope([],[],P_outer.A,P_outer.b),[1,2],'c');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       09-November-2022 
% Last update:   ---
% Last revision: 23-September-2024 (MW, integrate in polytope class)

% ------------------------------ BEGIN CODE -------------------------------

type = setDefaultValues({'outer'},varargin);

% convert the domain to an interval (outer approximation!)
if ~isa(dom,'interval')
    dom = interval(dom);
end

% normalize constraints
[A_,b_,Ae_,be_] = priv_normalizeConstraints(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,'A');

% compute intervals for normal vector and offset
A = interval(min(A_,[],1)',max(A_,[],1)');
b = interval(min(b_),max(b_));

% construct new halfspace containing all halfspaces in the list
A_reduced = center(A);
b_reduced = b - (A - A_reduced)'*dom;

if strcmp(type,'outer')
    P_out = polytope(A_reduced',supremum(b_reduced));
else
    P_out = polytope(A_reduced',infimum(b_reduced));
end

% ------------------------------ END OF CODE ------------------------------
