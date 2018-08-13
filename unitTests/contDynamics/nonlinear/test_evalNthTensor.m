function res = test_evalNthTensor(~)
% test_evalNthTensor - unit-test for the tensor generation and evaluation
%
% The solution from the evaluation of and generation with the toolbox
% funcitons "generateNthTensor" and "evalNthTensor" is compared to the
% evaluation using the corresponding closed-expression equation
%
% Syntax:  
%    res = test_evalNthTensor(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Niklas Kochdumper
% Written:      08-February-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------



% 2D-example --------------------------------------------------------------

% define test function
syms x y
f = sin(x)*cos(y+x)*exp(x*y);

% define test values
N = 10;
val = rand(2,10);

% compute derivatives
df_x = diff(f,x);
df_y = diff(f,y);

df_xx = diff(df_x,x);
df_xy = diff(df_x,y);
df_yy = diff(df_y,y);

df_xxx = diff(df_xx,x);
df_xxy = diff(df_xx,y);
df_xyy = diff(df_xy,y);
df_yyy = diff(df_yy,y);

df_xxxx = diff(df_xxx,x);
df_xxxy = diff(df_xxx,y);
df_xxyy = diff(df_xxy,y);
df_xyyy = diff(df_xyy,y);
df_yyyy = diff(df_yyy,y);

% evaluate function with formula
res_real = zeros(N,1);

for i = 1:N
   p = val(:,i);
   
   first = eval(subs(df_x,[x;y],p)) * p(1) + eval(subs(df_y,[x;y],p)) * p(2);
   second = 0.5 * eval(subs(df_xx,[x;y],p)) * p(1)^2 + ...
            eval(subs(df_xy,[x;y],p)) * p(1) * p(2) + ...
            0.5 * eval(subs(df_yy,[x;y],p)) * p(2)^2;
   third = 1/6 * eval(subs(df_xxx,[x;y],p)) * p(1)^3 + ...
           0.5 * eval(subs(df_xxy,[x;y],p)) * p(1)^2 * p(2) + ...
           0.5 * eval(subs(df_xyy,[x;y],p)) * p(1) * p(2)^2 + ...
           1/6 * eval(subs(df_yyy,[x;y],p)) * p(2)^3;
       
   fourth = 1/24 * eval(subs(df_xxxx,[x;y],p)) * p(1)^4 + ...
            1/6 * eval(subs(df_xxxy,[x;y],p)) * p(1)^3 * p(2) + ...
            1/4 * eval(subs(df_xxyy,[x;y],p)) * p(1)^2 * p(2)^2 + ...
            1/6 * eval(subs(df_xyyy,[x;y],p)) * p(1) * p(2)^3 + ...
            1/24 * eval(subs(df_yyyy,[x;y],p)) * p(2)^4;
       
   res_real(i) = first + second + third + fourth;
end

% evaluate function with method "evalNthTensor"
T = cell(3,1);
T{1} = generateNthTensor(f,[x;y],1);
T{2} = generateNthTensor(f,[x;y],2);
T{3} = generateNthTensor(f,[x;y],3);
T{4} = generateNthTensor(f,[x;y],4);

first = evalNthTensor(T{1},[x;y],1);
second = evalNthTensor(T{2},[x;y],2);
third = evalNthTensor(T{3},[x;y],3);
fourth = evalNthTensor(T{4},[x;y],4);

res_test = zeros(N,1);

for i = 1:N
    p = val(:,i);
    
    res_test(i) = eval(subs(first,[x;y],p)) + eval(subs(second,[x;y],p)) + ...
                  eval(subs(third,[x;y],p)) + eval(subs(fourth,[x;y],p));
end

% compare the results
res = 1;

for i = 1:N
    if abs(res_test(i)-res_real(i)) > 1e-14
       res = 0;
       break;
    end
end

%------------- END OF CODE --------------