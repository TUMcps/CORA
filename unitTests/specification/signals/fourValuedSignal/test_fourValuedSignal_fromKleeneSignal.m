function res = test_fourValuedSignal_fromKleeneSignal
% test_fourValuedSignal_fromKleeneSignal - unit test function of fromKleeneSignal
%
% Syntax:
%    res = test_fourValuedSignal_fromKleeneSignal
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
% Written:       20-February-2024
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
    sig = kleeneSignal(signals{i},signals{mod(i+1,length(signals))+1});
    expected = fourValuedSignal(sig,blank);
    actual = fourValuedSignal.fromKleeneSignal(sig);
    assertLoop(actual == expected,i)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
