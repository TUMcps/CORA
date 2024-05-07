function res = test_zonotope_ellipsoid
% test_zonotope_ellipsoid - unit test function of ellipsoid
%
% Syntax:
%    res = test_zonotope_ellipsoid
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% create zonotopes
Z1 = zonotope([-4, -3, -2; 1, 2, 3]);
Z2 = zonotope([1, 4, 2,1; -3, 2, -1,1]);

E1o = ellipsoid(Z1);
resvec(end+1) = true;
E1i = ellipsoid(Z1,'inner:norm');
resvec(end+1) = true;

Y1i = randPoint(E1i,2*dim(Z1),'extreme');

if ~contains(E1o,E1i) || ~all(contains(Z1,Y1i)) || ~all(contains(E1o,Y1i))
    resvec(end+1) = false;
end

E2o = ellipsoid(Z2,'outer:exact');
resvec(end+1) = true;
E2i = ellipsoid(Z2,'inner:exact');
resvec(end+1) = true;

Y2i = randPoint(E2i,2*dim(Z2),'extreme');
    
if ~contains(E2o,E2i) || ~all(contains(Z2,Y2i)) || ~all(contains(E2o,Y2i))
    resvec(end+1) = false;
end

% check point
c = [3;4;2];
Z = zonotope(c);
E = ellipsoid(Z);
[res,p] = representsa(E, 'point');
resvec(end+1) = res && all(p == c);

% gather results
res = all(resvec);
    
end


% ------------------------------ END OF CODE ------------------------------
