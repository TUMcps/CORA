function res = test_probZonotope_pyramid
% test_probZonotope_pyramid - unit test to check whether the probability of
%    intersection of a probabilistic zonotope with a polytope can be
%    appropriately computed by the pyramid function
%
% Syntax:
%    res = test_probZonotope_pyramid
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

% Authors:       Matthias Althoff
% Written:       17-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate probabilistic zonotope
Z1 = [0 1 0; 0 0 1];
Z2 = [0.6 1.2  ; 0.6 -1.2 ];
pZ3 = probZonotope(Z1,Z2,2);

% figure
% hold on
% plot(pZ3,'meshHide')

% polytope to be intersected
I = interval([-10; -5],[10; -3]);
P = polytope(I);

% provide mArray to determine accuracy of over-approximation
mArray = [3 2 1];

% compute intersection probability
intersectionProb = pyramid(pZ3,mArray,P); 

% saved result
intersectionProb_saved = 0.4867005538850111;

res = withinTol(intersectionProb,intersectionProb_saved,1e-8);

% ------------------------------ END OF CODE ------------------------------
