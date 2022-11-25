function res = testINTLAB_power_random(~)
% test_power_random - unit_test_function for comparing to IntLabV6
%
% Syntax:  
%    res = test_power_random (~)
%
% Inputs:
%    no 
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

% Author:       Dmitry Grebenyuk
% Written:      10-February-2016
% Last update:  
% Last revision:---

%------------- BEGIN CODE --------------
tol = 1e-2;
res = true;

try
    intvalinit('SharpIVmult');
catch
    res = false;
    disp('intvalinit failed');
    return;
end

format shortEng
format compact

a = 0;
b = 100;
min1 = (b-a).*rand(10000,1) + a;

a = 0;
b = 100;
delta1 = (b-a).*rand(10000,1) + a;

a = -5;
b = 5;
min2 = (b-a).*rand(10000,1) + a;

a = 0;
b = 5;
delta2 = (b-a).*rand(10000,1) + a;

int0 = interval(min1, min1 + delta1) .^ -3.3; %interval(min2, min2 + delta2);
int1 = infsup(min1, min1 + delta1) .^ -3.3; %infsup(min2, min2 + delta2);

i0 = infimum(int0);
i1 = inf(int1);
s0 = supremum(int0);
s1 = sup(int1);

[i0 , i1, i0 - i1];  %! Delete semicolumn to see [number, infimum in Cora, infimum in IntLab, diference]
[s0 , s1, s0 - s1];  %! Delete semicolumn to see [number, supremum in Cora, supremum in IntLab, diference]

bad_ones_min = find(abs(i0 - i1) > tol);
bad_ones_max = find(abs(s0 - s1) > tol);  

format long

if ( isempty(bad_ones_min) ~= true)
    disp('Infinums with diference > 0.000000001')
    disp('[number, infimum in Cora, infimum in IntLab, diference]')
    [bad_ones_min, i0(bad_ones_min), i1(bad_ones_min), i0(bad_ones_min) - i1(bad_ones_min)]
    disp(' ')
    res = false;
end

if ( isempty(bad_ones_max) ~= true)
    disp('Supremums with diference > 0.000000001')
    disp('LEGEND: [number, supremum in Cora, supremum in IntLab, diference]')
    [bad_ones_max, s0(bad_ones_max), s1(bad_ones_max), s0(bad_ones_max) - s1(bad_ones_max)]
    disp(' ')
    res = false;
end

disp('test_power_random successful');
disp(' ')

return;

