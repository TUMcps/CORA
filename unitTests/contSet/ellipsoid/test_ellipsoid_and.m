function res = test_ellipsoid_and
% test_ellipsoid_and - unit test function of and
%
% Syntax:
%    res = test_ellipsoid_and
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       26-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
% empty set

% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);
Ed1 = ellipsoid([ 4.2533342807136076 0.6346400221575308 ; 0.6346400221575309 0.0946946398147988 ], [ -2.4653656883489115 ; 0.2717868749873985 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);

% intersection with empty set
% res = (E_c{1}.E1 & E_e) == E_e;
    
% test ellipsoid (degenerate)
Eres = E1 & Ed1;
Yd = randPoint(Ed1,2);
for j=1:size(Yd,2)
    if contains(E1,Yd(:,j))
        assertLoop(contains(Eres,Yd(:,j)),j)
    end
end

% test ellipsoid (all zero)
Eres_0 = E1 & E0;
if contains(E1,E0.q)
    assert(rank(Eres_0)==0)
    assert(withinTol(Eres_0.q,E0.q,Eres_0.TOl))
else
    assert(representsa_(Eres_0,'emptySet',eps))
end

% test polytope
val = supportFunc(E1,unitvector(1,dim(E1)));
P = polytope(unitvector(1,dim(E1))',val);

Eres_h = E1 & P;

Y = randPoint(E1,2);
assert(all(contains(Eres_h,Y)))

% ------------------------------ END OF CODE ------------------------------
