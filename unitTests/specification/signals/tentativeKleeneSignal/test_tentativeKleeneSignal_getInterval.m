function res = test_tentativeKleeneSignal_getInterval
% test_tentativeKleeneSignal_getInterval - unit test function of getInterval
%
% Syntax:
%    res = test_tentativeKleeneSignal_getInterval
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       21-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
canBeTrueSig = pointSegmentSignal.indicator(stlInterval(1,3,true,false),true,false);
canBeFalseSig = pointSegmentSignal.indicator(stlInterval(2,4,false,true),true,false);
sig = tentativeKleeneSignal(canBeTrueSig,canBeFalseSig);

% test case definition
tt = fourValued.True;
uu = fourValued.Unknown;
ff = fourValued.False;
ii = fourValued.Inconclusive;
test_cases = {
    % {interval, expected}
    {stlInterval(1,4), fourValuedSignal.indicator(stlInterval(1,2),tt,ii).set(stlInterval(2,3,false),uu).set(stlInterval(3,4),ff)};
    {stlInterval(1,2), fourValuedSignal.indicator(stlInterval(1,2),tt,ii)};
    {stlInterval(2,3,false), fourValuedSignal.indicator(stlInterval(2,3,false),uu,ii)};
    {stlInterval(3,4), fourValuedSignal.indicator(stlInterval(3,4),ff,ii)};
    {stlInterval(1,2.5), fourValuedSignal.indicator(stlInterval(1,2),tt,ii).set(stlInterval(2,2.5,false,true),uu)};
    {stlInterval(0,inf), fourValuedSignal.indicator(stlInterval(1,2),tt,ii).set(stlInterval(2,3,false),uu).set(stlInterval(3,4),ff)};
    {stlInterval(), fourValuedSignal.uniformSignal(ii)};
};

% run tests
for i = 1:length(test_cases)
    interval = test_cases{i}{1};
    expected = test_cases{i}{2};
    actual = sig.getInterval(interval);
    assert(aux_intervalsEqual(expected,actual))
end

res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_intervalsEqual(sig1,sig2)
    fourValuedVals = enumeration('fourValued');
    for i = 1:length(fourValuedVals)
        int1 = sig1.findIntervals(fourValuedVals(i));
        int2 = sig2.findIntervals(fourValuedVals(i));
        assert(all(arrayfun(@(x,y) x == y,int1,int2)));
    end
    res = true;
end

% ------------------------------ END OF CODE ------------------------------
