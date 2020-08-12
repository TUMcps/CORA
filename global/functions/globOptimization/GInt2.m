function res = GInt2(var1, var2, funtion, eps)
% GInt2 - global optimization for intervals. The x, y have to be used as 
% variables in an inpunt function, but var1, var2 can be named what ever you
% wish.
%
% Syntax:  
%    res = GInt1(var1, var2, funtion)
%
% Inputs:
%    function - input funtion, string ' 1 +  x.^2 .* y.^2 + x '
%    var1, 2 - input variables, intervals.
%    eps - the precision, abs(res(i) - res(i - 1)) <= eps
%    
%
% Outputs:
%    res - an interval object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author:       Dmitry Grebenyuk
% Written:      07-July-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

%tic

% intitiate the var1, var2 as the x, y
x = var1;
y = var2;
% res_old is basicaly res(i-1) step. Evaluate the function using initial
% intervals
res_old = eval(funtion);
int = res_old;

% find out the most changing variable
x = split(x, [1, 1]);
resX = abs( rad(unite_ints(eval(funtion))) - rad(res_old) );
x = var1;
y = split(y, [1, 1]);
resY = abs( rad(unite_ints(eval(funtion))) - rad(res_old) );
y = var2;

% then subdivide that variable
if resX >= resY % x
    while 1

        % find the ordinal numbers of the minimum and the maximum of the
        % function
        k = find_ordinal_numbers(int);
        % devide the intervals containing the minimum and the maximum of the
        % function by two
        x = split(x, k);
        % evaluate the function intevalwise using new intervals
        int = eval(funtion);
        % unite the abtained result into one interval
        res = unite_ints(int);

        % when abs(res(i) - res(i - 1)) <= eps, it stops
        if (abs(infimum(res) - infimum(res_old)) <= eps) && (abs(supremum(res) - supremum(res_old)) <= eps)
            break
        else
            res_old = res;
        end

    end
    
else % y
    while 1

        k = find_ordinal_numbers(int);
        y = split(y, k);
        int = eval(funtion);

        res = unite_ints(int);

        if (abs(infimum(res) - infimum(res_old)) <= eps) && (abs(supremum(res) - supremum(res_old)) <= eps)
            break
        else
            res_old = res;
        end
    end
    
end

if 0
res_x = x;

x = var1;
y = var2;
res_old = eval(funtion);
int = res_old;
while 1

    k = find_ordinal_numbers(int);
    y = split(y, k);
    int = eval(funtion);

    res = unite_ints(int);

    if (abs(infimum(res) - infimum(res_old)) <= eps) && (abs(supremum(res) - supremum(res_old)) <= eps)
        break
    else
        res_old = res;
    end
end
res_y = y;

x_final_inf(k) = 0;
x_final_sup(k) = 0;
y_final_inf(k) = 0;
y_final_sup(k) = 0;

x_inf = infimum(res_x);
x_sup = supremum(res_x);
y_inf = infimum(res_y);
y_sup = supremum(res_y);

k = 1;
for j = 1:length(res_y)
    for i = 1:length(res_x)     
        %subs_i.type = '()';
        %subs_i.subs = {i};
        %subs_j.type = '()';
        %subs_j.subs = {j};
        x_final_inf(k) = x_inf(i);
        x_final_sup(k) = x_sup(i);
        y_final_inf(k) = y_inf(j);
        y_final_sup(k) =y_sup(j);
        %y_final(k) = res_y(j);
        %y1 = subsref(p1,subs);
        k = k + 1;
    end
end
x = interval(x_final_inf, x_final_sup);
y = interval(y_final_inf, y_final_sup);

int = eval(funtion);
min_int = infimum(int);
max_int = supremum(int);
res = interval(min(min_int), max(max_int));
end

%global2_t3 = toc

%------------- END OF CODE --------------
end