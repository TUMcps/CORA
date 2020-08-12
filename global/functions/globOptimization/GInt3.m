function res = GInt3(var1, var2, var3, funtion, eps)
% GInt3 - global optimization for intervals. The x, y, z have to be used as
% variables in an inpunt function, but var1, var2, var3 can be named what ever you
% wish.
%
% Syntax:  
%    res = GInt1(var1, var2, var3, funtion)
%
% Inputs:
%    function - input funtion, string ' 1 +  x.^2 .* y.^2 .* z.^3 + x '
%    var1, 2, 3 - input variables, intervals.
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
% Written:      12-July-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

%tic

% intitiate the var1, var2 as the x, y, z
x = var1;
y = var2;
z = var3;
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

z = split(z, [1, 1]);
resZ = abs( rad(unite_ints(eval(funtion))) - rad(res_old) );
z = var3;

% then subdivide that variable
if (resX >= resY) && (resX >= resZ)     % x
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
    
elseif resY >= resX && resY >= resZ     % y
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
    
else    % z
    while 1

        k = find_ordinal_numbers(int);
        z = split(z, k);
        int = eval(funtion);

        res = unite_ints(int);

        if (abs(infimum(res) - infimum(res_old)) <= eps) && (abs(supremum(res) - supremum(res_old)) <= eps)
            break
        else
            res_old = res;
        end
    end
end

%global2_t3 = toc

%------------- END OF CODE --------------
end
