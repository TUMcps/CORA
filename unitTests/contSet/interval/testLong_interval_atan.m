function res = testLong_interval_atan
% testLong_interval_atan - unit test function of sine for intervals,
%    overloaded 'atan()' function for intervals
%
% Syntax:
%    res = testLong_interval_atan
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
% See also: mtimes

% Authors:       Mark Wetzlinger
% Written:       08-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Random tests
numberRandTests = 10000;

% all values must be between -pi/2 and pi/2
res_rand(1) = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    b = atan(a);
    if b.inf < -pi/2 || b.sup > pi/2
        res_rand(1) = false;
        break
    end
end

% all x<0 map to atan(x)<0, all x>0 map to atan(x)>0
res_rand(2) = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    b = atan(a);
    if b.inf >= 0 || b.sup <= 0
        res_rand(2) = false;
        break
    end
end

% final test result
res = all(res_rand);

% ------------------------------ END OF CODE ------------------------------
