function int = plus(summand1,summand2)
% plus - Overloaded '+' operator for STL intervals
%
% Syntax:
%    int = plus(summand1,summand2)
%
% Inputs:
%    summand1 - stlInterval object
%    summand2 - stlInterval object
%
% Outputs:
%    int - interval
%
% Example:
%    summand1 = stlInterval(1,2);
%    summand2 = stlInterval(0,1,false,false);
%    summand1 + summand2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine the interval object
[res,summand] = findClassArg(summand1,summand2,'stlInterval');

% check empty set
if representsa(res,'emptySet') || representsa(summand,'emptySet')
    int = stlInterval();
    return;
end

% stl + stl
if isa(summand,'stlInterval')
    lb = res.lower + summand.lower;
    ub = res.upper + summand.upper;
    lc = res.leftClosed && summand.leftClosed;
    rc = res.rightClosed && summand.rightClosed;
    int = stlInterval(lb,ub,lc,rc);

% stl + interval
elseif isa(summand,'interval')
    if dim(summand) ~= 1
        throw(CORAerror('CORA:dimensionMismatch',res,summand));
    end
    lb = res.lower + infimum(summand);
    ub = res.upper + supremum(summand);
    int = stlInterval(lb,ub,res.leftClosed,res.rightClosed);

% stl + numeric
elseif isa(summand,'numeric')
    lb = res.lower + summand;
    ub = res.upper + summand;
    int = stlInterval(lb,ub,res.leftClosed,res.rightClosed);
else
    throw(CORAerror('CORA:noops',summand1,summand2));
end

% ------------------------------ END OF CODE ------------------------------
