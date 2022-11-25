function res = test_interval_enclosePoints
% test_interval_enclosePoints - unit test function of
%    enclosePoints
%
% Syntax:  
%    res = test_interval_enclosePoints
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

% Author:       Mark Wetzlinger
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% point cloud
pts = [1 4 2 5 3 2 4 3 2 5; ...
       6 8 9 7 6 9 8 6 8 7];

% init enclosing interval (has to be smaller than -1 to 1)
I = interval.enclosePoints(pts);

% solution
I_true = interval([1; 6],[5; 9]);

% check with correct solution
res = I == I_true;


if res
    disp('test_enclosePoints successful');
else
    disp('test_enclosePoints failed');
end

%------------- END OF CODE --------------
