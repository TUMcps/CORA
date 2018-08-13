function res = minus(minuend,subtrahend)
% plus - Overloaded '-' operator for intervals
%
% Syntax:  
%    res = plus(summand1,summand2)
%
% Inputs:
%    minuend - interval or numerical value
%    subtrahend - interval or numerical value
%
% Outputs:
%    res - interval
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      25-June-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%init 
res = interval();

%Find an interval object
%Is minuend an interval?
if isa(minuend,'interval')
    %Is subtrahend an interval?
    if isa(subtrahend,'interval')
        %Calculate infimum and supremum
        res.inf = minuend.inf - subtrahend.sup;
        res.sup = minuend.sup - subtrahend.inf;
    else
        %Calculate infimum and supremum
        res.inf = minuend.inf - subtrahend;
        res.sup = minuend.sup - subtrahend;       
    end
else
    %minuend must be a particular value
    %Calculate infimum and supremum
    res.inf = minuend - subtrahend.sup;
    res.sup = minuend - subtrahend.inf;
end

%------------- END OF CODE --------------