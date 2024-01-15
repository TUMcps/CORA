function res = test_polytope_enclosePoints
% test_polytope_enclosePoints - unit test function of enclosePoints
%
% Syntax:
%    res = test_polytope_enclosePoints
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       04-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D case
p = [-2 -4 1 4 3];
P = polytope.enclosePoints(p);
% check that all points are contained
res = all(contains(P,p));


% 2D, unit square
p = [1 0; 0 1; -1 0; 0 -1]';
P = polytope.enclosePoints(p);
% check that all points are contained
res(end+1,1) = all(contains(P,p));

% 2D, add a point that does not matter
p = [1 0; 0 1; -1 0; 0 -1; 0 0]';
P_ = polytope.enclosePoints(p);
% should be equal to previous polytope
res(end+1,1) = P == P_;

% 2D, degenerate case
p = [1 0; 0 0; -1 0]';
P = polytope.enclosePoints(p);
res(end+1,1) = all(contains(P,p));


% 4D, points in the unit box
p = [0.3 -0.9  0.7  0.1 ;
     0.2 -0.3  0.4  0.7 ;
    -0.5  0.7 -0.8 -0.4 ;
     0.8  0.2 -0.3  0.9 ;
    -0.9  0.1  0.5 -0.8 ;
    -0.1 -0.6  0.9  0.9 ]';
P = polytope.enclosePoints(p);
% ...has to be contained in the interval
I = interval(-ones(4,1),ones(4,1));
res(end+1,1) = contains(I,P);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
