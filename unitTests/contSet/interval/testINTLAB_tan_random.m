function res = testINTLAB_tan_random()
% testINTLAB_tan_random - unit_test_function for comparing to IntLabV6
%
% Syntax:
%    res = testINTLAB_tan_random
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

% Authors:       Dmitry Grebenyuk, Matthias Althoff
% Written:       05-February-2016
% Last update:   09-July-2021 (MA, improved the unit test and added documentation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set numerical tolerance
tol = 1e-8;

% initialize result
res = true;

% initialize sharp computation of intervals
try
    intvalinit('SharpIVmult');
catch
    res = false;
    disp('intvalinit failed');
    return;
end

% create left limit
a = -4*pi;
b = 4*pi;
min = (b-a).*rand(10000,1) + a;

% create right limit
a = 0;
b = 3*pi;
delta = (b-a).*rand(10000,1) + a;

% instantiate intervals of CORA and INTLAB and compute results
int0 = tan(interval(min, min + delta));
int1 = tan(infsup(min, min + delta));

% return infimum and supremum using CORA and INTLAB
i0 = infimum(int0);
i1 = inf(int1);
s0 = supremum(int0);
s1 = sup(int1);

% find indices for which the result differs more than the accepted
% tolerance
bad_ones_min = find(abs(i0 - i1) > tol);
bad_ones_max = find(abs(s0 - s1) > tol);

format long

% display incorrect results for the infimum
if ( isempty(bad_ones_min) ~= true)
    disp('Infinums with difference > 0.000000001')
    disp('[number, infimum in Cora, infimum in IntLab, diference]')
    [bad_ones_min, i0(bad_ones_min), i1(bad_ones_min), i0(bad_ones_min) - i1(bad_ones_min)]
    disp(' ')
    res = false;
end

% display incorrect results for the supremum
if ( isempty(bad_ones_max) ~= true)
    disp('Supremums with difference > 0.000000001')
    disp('LEGEND: [number, supremum in Cora, supremum in IntLab, diference]')
    [bad_ones_max, s0(bad_ones_max), s1(bad_ones_max), s0(bad_ones_max) - s1(bad_ones_max)]
    disp(' ')
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
