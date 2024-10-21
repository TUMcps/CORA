function res = test_kleeneSignal_not
% test_kleeneSignal_not - unit test function of not
%
% Syntax:
%    res = test_kleeneSignal_not
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
% Written:       19-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
tt = true;
ff = false;
signals = {
    kleeneSignal(pointSegmentSignal([0,1,2,2.5],[tt,tt,ff,tt,tt,ff,ff,tt]),pointSegmentSignal.indicator(stlInterval(1.5,2.25,tt,ff),tt,ff));
    kleeneSignal(pointSegmentSignal([0,1,3,4,9.5],[tt,ff,ff,tt,ff,ff,tt,tt,tt,ff]),pointSegmentSignal.indicator(stlInterval(5,10),tt,ff));
};
kleeneVals = [kleene.True,kleene.Unknown,kleene.False];

for i = 1:length(signals)
    sig = signals{i};
    actual = ~sig;
    sigInt = arrayfun(@(match) sig.findIntervals(match),kleeneVals,'UniformOutput',false);
    actualInt = arrayfun(@(match) actual.findIntervals(match),kleeneVals,'UniformOutput',false);
    assertLoop(compareIntervals(sigInt{1},actualInt{3}),i) % true becomes false
    assertLoop(compareIntervals(sigInt{2},actualInt{2}),i) % unknown stays unknown
    assertLoop(compareIntervals(sigInt{3},actualInt{1}),i) % false becomes true
end

res = true;

% ------------------------------ END OF CODE ------------------------------
