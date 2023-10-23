function res = test_zonotope_enlarge
% test_zonotope_enlarge - unit test function of enlarge
%
% Syntax:
%    res = test_zonotope_enlarge
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
% Written:       26-July-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% create zonotope
c = [-4; 1];
G = [-3, -2, -1; 2, 3, 4];
Z = zonotope(c,G);

% compute enlarged zonotope
Z_ = enlarge(Z,[2;1.5]);

% obtain center and generator matrix
c_ = Z_.c;
G_ = Z_.G;

% true result
true_c = [-4; 1];
true_G = [-6, -4, -2; ...
            3, 4.5, 6];

% check result
resvec(end+1) = compareMatrices(c_,true_c) && compareMatrices(G_,true_G);

% check with scalar factor
Z_ = enlarge(Z,-3);

% obtain center and generator matrix
c_ = Z_.c;
G_ = Z_.G;

% true result
true_c = [-4; 1];
true_G = [9, 6, 3; ...
            -6, -9, -12];

% check result
resvec(end+1) = compareMatrices(c_,true_c) && compareMatrices(G_,true_G);

% check empty generator matrix
Z = zonotope([1;2]);
Z = enlarge(Z,2);
resvec(end+1) = isempty(Z.G);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
