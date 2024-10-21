function res = test_kleeneSignal_fromBoolSignal
% test_kleeneSignal_fromBoolSignal - unit test function of fromBoolSignal
%
% Syntax:
%    res = test_kleeneSignal_fromBoolSignal
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
    pointSegmentSignal(0,[tt,ff]);
    pointSegmentSignal(0,[ff,tt]);
    pointSegmentSignal([0,1,2],[ff,tt,tt,ff,tt,tt]);
    pointSegmentSignal([0,1,2],[ff,tt,ff,tt,ff,tt]);
    pointSegmentSignal([0,2,3],[tt,tt,ff,ff,ff,tt]);
};
blank = pointSegmentSignal(0,[ff ff]);

% run tests
for i = 1:length(signals)
    sig = signals{i};
    expected = kleeneSignal(sig,blank);
    actual = kleeneSignal.fromBoolSignal(sig);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
