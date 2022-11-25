function res = testLongDuration_conZonotope_interval
% testLongDuration_conZonotope_interval - unit test function for the
%    calculation of a bounding box of a constrained zonotope object
%
% Syntax:  
%    res = testLongDuration_conZonotope_interval
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
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      22-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;

% TEST 1: Random Test 3D --------------------------------------------------

% Generate random conZonotope object
points = rand(3,100);
ind = convhulln(points');
ind = unique(ind(:,1),'stable');
V = points(:,ind);

poly = mptPolytope(V');
cZono = conZonotope(poly);

% calculate interval
int = interval(cZono);

% compare with ground-truth for the vertices
V = vertices(cZono);
int_ = interval(min(V,[],2),max(V,[],2));

for i = 1:length(int)
   if abs(infimum(int_(i)) - infimum(int(i))) > 1e-10 || abs(supremum(int_(i)) - supremum(int(i))) > 1e-10
      file_name = strcat('testLongDuration_conZonotope_interval_1_', ...
                         datestr(now,'mm-dd-yyyy_HH-MM'));
                  
      file_path = fullfile(coraroot(), 'unitTests', 'failedTests', ...
                           file_name);
                           
      save(file_path, 'cZono')
      error('Test 1 failed!'); 
   end
end


res = true;


%------------- END OF CODE --------------