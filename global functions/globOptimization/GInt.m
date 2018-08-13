function res = GInt(var1, var2, funtion)
% GOP - global optimization
%
% Syntax:  
%    res = GOP(funtion)
%
% Inputs:
%    function - input funtion, char
%    var1, 2, 3 - input variables, intervals
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


clc
eps = 0.1;

x = var1;
y = var2;
res_old = eval(funtion);
int = res_old;

while 1

    k = find_ordinal_numbers(int);
    x = split(x, k);
    int = eval(funtion);

    min_int = infimum(int);
    max_int = supremum(int);
    res = interval(min(min_int), max(max_int));

    if (abs(infimum(res) - infimum(res_old)) <= eps) && (abs(supremum(res) - supremum(res_old)) <= eps)
        break
    else
        res_old = res;
    end

end
res_x = x;

x = var1;
y = var2;
res_old = eval(funtion);
int = res_old;
while 1

    k = find_ordinal_numbers(int);
    y = split(y, k);
    int = eval(funtion);

    min_int = infimum(int);
    max_int = supremum(int);
    res = interval(min(min_int), max(max_int));

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

%------------- END OF CODE --------------
end