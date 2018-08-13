function res = speed_error_test( ~ )
% speed_error_test - a universal unit test function for measuring
% implementation time and comparing errors with IntLab.
%
% Syntax:  
%    res = speed_error_test (~)
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
% Written:      24-February-2016
% Last update:  
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

intvalinit('SharpIVmult');
format long

% the itervation of the test
iterations = 1e+1;

for i = 1:iterations

    % the range of min1
    a = -100;
    b = 100;
    min1 = (b-a).*rand(100,100) + a;
    
    % the range of delta1
    a = 0;
    b = 100;
    delta1 = (b-a).*rand(100,100) + a;
    
    % the range of min2
    a = -10;
    b = 10;
    min2 = (b-a).*rand(100,100) + a;

    % the range of delta2
    a = 0;
    b = 10;
    delta2 = (b-a).*rand(100, 100) + a;


    int0pre = interval(min1, min1 + delta1);    % [min1, min1 + delta1] in CORA's interval
    int1pre = infsup(min1, min1 + delta1);      % [min1, min1 + delta1] in IntLab's interval
    int0pre2 = interval(min2, min2 + delta2);   % [min2, min2 + delta2] in CORA's interval
    int1pre2 = infsup(min2, min2 + delta2);     % [min2, min2 + delta2] in IntLab's interval
    
    % calculation of the speed of implementation 
    tic
    int0 = tan(int0pre);% .^ int0pre2;   % insert a test function here. For example int0 = sin(int0pre) or int0 = int0pre .* int0pre2
    t_CORA(i) = toc;
    tic
    int1 = tan(int1pre);% .^ int1pre2;   % insert a test function here too. For example int1 = sin(int1pre) or int1 = int1pre .* int1pre2
    t_INTLAB(i) = toc;
    
    
    % calculation of the absolute and relative errors
    i0 = infimum(int0);
    i1 = inf(int1);
    s0 = supremum(int0);
    s1 = sup(int1);
    dev_den = (s0 - i0);
    dev_num = max(abs(i0 - i1), abs(s0 - s1));
    [i0 , i1, (i0 - i1)];  %! Delete semicolumn to see [number, infimum in Cora, infimum in IntLab, diference]
    [s0 , s1, s0 - s1];
    
    error_rel(i) = max(dev_num(:) ./ dev_den(:));
    error_abs(i) = max(dev_num(:));

    %tic
    %cosh(int0);
    %t_CORA(i) = toc;
    %tic
    %cosh(int1);
    %t_IntLab(i) = toc;

end
clc
disp('tan')         % the name of testing function here
t_CORA;
time_CORA = sum(t_CORA)/iterations
t_INTLAB;
time_INTLAB = sum(t_INTLAB)/iterations
error_rel;
max_relative_error = max(error_rel)
error_abs;
max_absolute_error = max(error_abs)

tol = 1e-9;

bad_ones_min = find(abs(i0 - i1) > tol);
bad_ones_max = find(abs(s0 - s1) > tol);  

format long

if ( isempty(bad_ones_min) ~= true)
    disp('Infinums with difference > 0.000000001')
    disp('[number, infimum in Cora, infimum in IntLab, difference]')
    [bad_ones_min, i0(bad_ones_min), i1(bad_ones_min), i0(bad_ones_min) - i1(bad_ones_min)]
    disp(' ')
    res = false;
end

if ( isempty(bad_ones_max) ~= true)
    disp('Supremums with difference > 0.000000001')
    disp('LEGEND: [number, supremum in Cora, supremum in IntLab, difference]')
    [bad_ones_max, s0(bad_ones_max), s1(bad_ones_max), s0(bad_ones_max) - s1(bad_ones_max)]
    disp(' ')
    res = false;
end

disp('speed_error_test is finished');
disp(' ')

end

%------------- END OF CODE --------------
