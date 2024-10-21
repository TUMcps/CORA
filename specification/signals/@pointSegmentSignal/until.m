function sig = until(lhs,interval,rhs)
% until - compute the until combination of two signals
%
% Syntax:
%    sig = until(lhs,interval,rhs)
%
% Inputs:
%    lhs - left-hand side signal
%    interval - time interval of the until operator
%    rhs - right-hand side signal
%
% Outputs:
%    sig - the combined signal
%
% Example:
%    lhs = pointSegmentSignal(0,[true true]);
%    rhs = pointSegmentSignal([0 4 7],[false false true true true false]);
%    int = stlInterval(1, 2);
%    sig = until(lhs,int,rhs);
%    assert(sig == pointSegmentSignal([0 2 6],[false false true true true false]));
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       12-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find the positive intervals
lhsInt = lhs.findIntervals();
rhsInt = rhs.findIntervals();

% set to true wherever we find a witness for until
sig = aux_intervalUntil(lhsInt,interval,rhsInt);

end


% Auxiliary functions -----------------------------------------------------

function out = aux_intervalUntil(lhsInt,untilInt,rhsInt)
% compute the until combination of two signals given their positive intervals

% set output to false by default
out = pointSegmentSignal.indicator(stlInterval(0,inf),false,false);

if untilInt.contains(0)
    [ub,rc] = supremum(untilInt);
    untilIntNoZero = stlInterval(0,ub,false,rc);
else
    untilIntNoZero = untilInt;
end
for i=1:length(lhsInt)
    for j=1:length(rhsInt)
        posShift = ((rhsInt(j) & rightClosure(lhsInt(i))) - untilIntNoZero) & leftClosure(lhsInt(i));

        if ~representsa(posShift,'emptySet')
            out = out.set(posShift,true);
        end
    end
end

% if 0 delay is allowed, then until is true where rhs is true
% regardless of lhs (since we use the strict semantics)
if untilInt.contains(0)
    for j=1:length(rhsInt)
        out = out.set(rhsInt(j),true);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
