function res = testINTLAB_sin_random(~)
% test_sin_random - unit_test_function for comparing to IntLabV6
%
% Syntax:  
%    res = test_sin_random (~)
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
% Written:      05-February-2016
% Last update:  
% Last revision:---

%------------- BEGIN CODE --------------
tol = 1e-9;
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

a = -4*pi;
b = 4*pi;
min = (b-a).*rand(10000,1) + a;

a = 0;
b = 3*pi;
delta = (b-a).*rand(10000,1) + a;

angl0 = min;
angl1 = min + delta;

int0 = sin(interval(min, min + delta));
int1 = sin(infsup(min, min + delta));

i0 = infimum(int0);
i1 = inf(int1);
s0 = supremum(int0);
s1 = sup(int1);

%diference_min = infimum(int0) - inf(int1);
%diference_max = supremum(int0) - sup(int1);

[i0 , i1, i0 - i1];  %! Delete semicolumn to see [number, infimum in Cora, infimum in IntLab, diference]
[s0 , s1, s0 - s1];  %! Delete semicolumn to see [number, supremum in Cora, supremum in IntLab, diference]

bad_ones_min = find(abs(i0 - i1) > 0.000000001);
bad_ones_max = find(abs(s0 - s1) > 0.000000001);

format long

%disp('PROBLEM ONES - > 0.000000001 diference')

%for i = 1:10000.
%    if (i0(i) - i1(i) > 0.000000001 ) && ( i0(i) - i1(i) ~= Inf )
%        [i, i0(i), i1(i), i0(i) - i1(i)]
%        res = 0;
%    end
%end

%for i = 1:10000.
%    if (s0(i) - s1(i) > 0.000000001 ) && ( s0(i) - s1(i) ~= Inf )
%        [i, s0(i), s1(i), s0(i) - s1(i)]
%        res = 0;
%    end
%end

if ( isempty(bad_ones_min) ~= true)
    disp('Infinums with diference > 0.000000001')
    disp('[number, infimum in Cora, infimum in IntLab, diference]')
    [bad_ones_min, i0(bad_ones_min), i1(bad_ones_min), i0(bad_ones_min) - i1(bad_ones_min), 180*angl0(bad_ones_min)/pi, 180*angl1(bad_ones_min)/pi]
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

disp('test_sin_random successful');
disp(' ')

return;
