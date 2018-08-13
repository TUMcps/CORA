function res = GInt1(var1, funtion, eps)
% GInt1 - global optimization for intervals. The x has to be used as a
% variable in an inpunt function, but var1 can be named what ever you
% wish.
%
% Syntax:  
%    res = GInt1(var1, funtion)
%
% Inputs:
%    function - input funtion, string ' 1 +  x.^2 '
%    var1 - an input variable, interval
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

% intitiate the var1 as the x
x = var1;

% res_old is basicaly res(i-1) step. Evaluate the function using initial
% intervals
res_old = eval(funtion);
int = res_old;

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

% measure time
%global2_t3 = toc

%------------- END OF CODE --------------
end

