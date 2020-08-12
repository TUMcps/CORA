function res = test_in
% test_in - unit test function of in
%
% Syntax:  
%    res = test_in
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
% Written:      27-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

I = interval([-3;-2],[5;4]);
Z_in = zonotope([0.5, 2, 1;
                 0,   1,-0.7]);
Z_out = zonotope([6.5, 2, 1;
                 -3,   1,-0.7]);

res_in = in(I, Z_in);
res_out = in(I, Z_out); % should be out, but interval overapprox is in

res = res_in && ~res_out;


if res
    disp('test_in successful');
else
    disp('test_in failed');
end

%------------- END OF CODE --------------