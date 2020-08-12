function res = test_zonotope_volume
% test_zonotope_volume - unit test function of volume
%
% Syntax:  
%    res = test_zonotope_volume
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  01-May-2020 (MW, add second case)
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% obtain result
vol = volume(Z1);

% true result 1
true_vol = 80;

res_true = vol == true_vol;

% compare to interval
I1 = interval(Z1);
volInt = volume(I1);
% convert back to zonotope
IZ1 = zonotope(I1);
volIntZon = volume(IZ1); % has to be equal to interval volume

res_int1 = vol < volInt;
res_int2 = volIntZon == volInt;

% check results
res = res_true && res_int1 && res_int2;


if res
    disp('test_zonotope_volume successful');
else
    disp('test_zonotope_volume failed');
end

%------------- END OF CODE --------------
