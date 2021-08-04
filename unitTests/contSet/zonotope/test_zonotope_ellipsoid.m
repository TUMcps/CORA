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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;

% create zonotopes
Z1 = zonotope([-4, -3, -2; 1, 2, 3]);
Z2 = zonotope([1, 4, 2,1; -3, 2, -1,1]);

E1o = ellipsoid(Z1);
E1i = ellipsoid(Z1,'i:norm');

Y1i = randPoint(E1i,2*dim(Z1),'extreme');

if ~in(E1o,E1i) || ~in(Z1,Y1i) || ~in(E1o,Y1i)
    res = false;
end

E2o = ellipsoid(Z2,'o:exact');
E2i = ellipsoid(Z2,'i:exact');

Y2i = randPoint(E2i,2*dim(Z2),'extreme');
    
if ~in(E2o,E2i) || ~in(Z2,Y2i) || ~in(E2o,Y2i)
    res = false;
end

if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end

%------------- END OF CODE --------------
