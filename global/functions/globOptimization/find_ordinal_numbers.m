function res = find_ordinal_numbers(input)
% find_ordinal_numbers - finds the minimum and maximum of ordinal numbers of an interval
%
% Syntax:  
%    res = find_ordinal_numbers(input)
%
% Inputs:
%  
%   input - an input interval
%
% Outputs:
%    res - minimum and maximum of ordinal numbers
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk
% Written:      06-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

inf = input.inf;
sup = input.sup;

min_int = inf(1);
k_min = 1;
max_int = sup(1);
k_max = 1;

for i = 2: length(inf)
    if min_int > inf(i)
        k_min = i;
        min_int = inf(i);
    end
    if max_int < sup(i)
        k_max = i;
        max_int = sup(i);
    end
end

res = [k_min, k_max];

%------------- END OF CODE --------------
end

