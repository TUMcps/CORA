function res = test_mptPolytope_project
% test_mptPolytope_project - unit test function for projection of polytopes
%
% Syntax:  
%    res = test_mptPolytope_project()
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

% Author:       Niklas Kochdumper
% Written:      21-December-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

% Analytical Tests --------------------------------------------------------

% define polytope
A = [-1 0; 2 4; 1 -2];
b = [-1; 14; -1];

poly = mptPolytope(A,b);

% project to dimension 1
poly_ = project(poly,1);
V = vertices(poly_);

if ~all(ismember(V,[1 3]))
   error('Analytical test for "mptPolytope/project" failed!');
end

% project to dimension 2
poly_ = project(poly,2);
V = vertices(poly_);

if ~all(ismember(V,[1 3]))
   error('Analytical test for "mptPolytope/project" failed!');
end



res = true;
    