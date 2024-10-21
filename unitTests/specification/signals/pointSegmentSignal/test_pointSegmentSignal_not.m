function res = test_pointSegmentSignal_not
% test_pointSegmentSignal_not - unit test function of not
%
% Syntax:
%    res = test_pointSegmentSignal_not
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
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% signals
tt = true;
ff = false;
signals = {
    pointSegmentSignal([0,1,2,2.5], [tt,tt,ff,tt,tt,ff,ff,tt]);
    pointSegmentSignal([0,1,2,2.5], [tt,ff,ff,tt,ff,ff,ff,tt]);
    pointSegmentSignal(0, [ff,ff]);
    pointSegmentSignal(0, [tt,ff]);
    pointSegmentSignal([0,0.5,1,1.5,2], [tt,ff,tt,ff,tt,ff,tt,ff,tt,tt]);
    pointSegmentSignal([0,0.75,1,1.25,2], [tt,tt,ff,ff,tt,ff,ff,tt,tt,ff]);    
};

% run tests
for i = 1:length(signals)
    sig = signals{i};
    expected = pointSegmentSignal(sig.timePoints,~sig.values);
    actual = ~sig;
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
