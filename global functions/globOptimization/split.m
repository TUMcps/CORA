function res = split(input, number)
% split - splits an interval into two intervals
%
% Syntax:  
%    res = split(input, number)
%
% Inputs:
%  
%   input - an input interval
%   number - is an ordinal number of the interval subset
%
% Outputs:
%    res - an interval
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

if number(1) == number(2)
    new_inf = [inf(1:number(1)), (sup(number(1)) + inf(number(1)))/2, inf(number(1)+1:end) ];   
    new_sup = [sup(1:number(1)-1), (sup(number(1)) + inf(number(1)))/2, sup(number(1):end) ];
    
elseif number(1) < number(2)
    new_inf = [inf(1:number(1)), (sup(number(1)) + inf(number(1)))/2, inf(number(1)+1:number(2)), (sup(number(2)) + inf(number(2)))/2, inf(number(2)+1:end) ];   
    new_sup = [sup(1:number(1)-1), (sup(number(1)) + inf(number(1)))/2, sup(number(1):number(2)-1), (sup(number(2)) + inf(number(2)))/2, sup(number(2):end) ];
    
else
    new_inf = [inf(1:number(2)), (sup(number(2)) + inf(number(2)))/2, inf(number(2)+1:number(1)), (sup(number(1)) + inf(number(1)))/2, inf(number(1)+1:end) ];   
    new_sup = [sup(1:number(2)-1), (sup(number(2)) + inf(number(2)))/2, sup(number(2):number(1)-1), (sup(number(1)) + inf(number(1)))/2, sup(number(1):end) ];
end

res = interval(new_inf, new_sup);

%------------- END OF CODE --------------

end

