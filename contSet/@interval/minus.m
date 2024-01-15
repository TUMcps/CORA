function res = minus(minuend,subtrahend)
% minus - Overloaded '-' operator for intervals
%
% Syntax:
%    res = minus(minuend,subtrahend)
%
% Inputs:
%    minuend - interval or numerical value
%    subtrahend - interval or numerical value
%
% Outputs:
%    res - interval
%
% Example:
%    minuend = interval([-2;1],[3;2]);
%    subtrahend = interval([-0.5;0.2],[0.2;0.6]);
%    minuend - subtrahend
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       25-June-2015
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%Find an interval object
%Is minuend an interval?
if isa(minuend,'interval')
    %init 
    res = minuend;
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
    %init 
    res = subtrahend;
    %minuend must be a particular value
    %Calculate infimum and supremum
    res.inf = minuend - subtrahend.sup;
    res.sup = minuend - subtrahend.inf;
end

% ------------------------------ END OF CODE ------------------------------
