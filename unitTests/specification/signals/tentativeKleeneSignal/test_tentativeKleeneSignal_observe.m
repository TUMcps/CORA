function res = test_tentativeKleeneSignal_observe
% test_tentativeKleeneSignal_observe - unit test function of observe
%
% Syntax:
%    res = test_tentativeKleeneSignal_observe
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
sig = tentativeKleeneSignal.emptySignal();
expectedT = pointSegmentSignal.indicator(stlInterval(),false,false);
expectedF = pointSegmentSignal.indicator(stlInterval(),false,false);

tt = kleene.True;
uu = kleene.Unknown;
ff = kleene.False;
observations = {
    % {interval, value}
    {stlInterval(0,1), tt};
    {stlInterval(1,2), uu};
    {stlInterval(2,3), ff};
    {stlInterval(3,4), tt};
    {stlInterval(0,3,false), ff};
    {stlInterval(0), tt};
    {stlInterval(2.5,3.5,false,true), uu};
};

for i = 1:length(observations)
    interval = observations{i}{1};
    value = observations{i}{2};
    sig = sig.observe(interval,value);
    switch value
        case kleene.True
            expectedT = expectedT.set(interval,true);
        case kleene.Unknown
            expectedT = expectedT.set(interval,true);
            expectedF = expectedF.set(interval,true);
        case kleene.False
            expectedF = expectedF.set(interval,true);
    end
    assertLoop(sig.canBeTrue == expectedT,i)
    assertLoop(sig.canBeFalse == expectedF,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
