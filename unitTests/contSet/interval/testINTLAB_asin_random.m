function res = testINTLAB_asin_random(~)
% testINTLAB_asin_random - unit_test_function for comparing to IntLabV6
%
% Syntax:
%    res = testINTLAB_asin_random
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

% Authors:       Dmitry Grebenyuk
% Written:       06-February-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true;

try
    intvalinit('SharpIVmult');
catch
    assert(false);
    disp('intvalinit failed');
    return;
end

format shortEng
format compact

% start test
a = -1;
b = 0;
min = (b-a).*rand(10000,1) + a;

a = 0;
b = 1;
delta = (b-a).*rand(10000,1) + a;

int0 = asin(interval(min, min + delta));
int1 = asin(infsup(min, min + delta));

i0 = infimum(int0);
i1 = inf(int1);
s0 = supremum(int0);
s1 = sup(int1);

[i0 , i1, i0 - i1];  %! Delete semicolumn to see [number, infimum in Cora, infimum in IntLab, diference]
[s0 , s1, s0 - s1];  %! Delete semicolumn to see [number, supremum in Cora, supremum in IntLab, diference]

bad_ones_min = find(abs(i0 - i1) > 0.00000001);
bad_ones_max = find(abs(s0 - s1) > 0.00000001);  

format long

if ( isempty(bad_ones_min) ~= true)
    disp('Infinums with diference > 0.000000001')
    disp('[number, infimum in Cora, infimum in IntLab, diference]')
    [bad_ones_min, i0(bad_ones_min), i1(bad_ones_min), i0(bad_ones_min) - i1(bad_ones_min)]
    disp(' ')
    assert(false);
end

if ( isempty(bad_ones_max) ~= true)
    disp('Supremums with diference > 0.000000001')
    disp('LEGEND: [number, supremum in Cora, supremum in IntLab, diference]')
    [bad_ones_max, s0(bad_ones_max), s1(bad_ones_max), s0(bad_ones_max) - s1(bad_ones_max)]
    disp(' ')
    assert(false);
end

% ------------------------------ END OF CODE ------------------------------
