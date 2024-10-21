function res = test_tentativeKleeneSignal_setInvalid
% test_tentativeKleeneSignal_setInvalid - unit test function of setInvalid
%
% Syntax:
%    res = test_tentativeKleeneSignal_setInvalid
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

sig = sig.setInvalid(stlInterval(2,3,false));
assert(sig.canBeTrue == pointSegmentSignal.indicator(stlInterval(1,2),true,false));
assert(sig.canBeFalse == pointSegmentSignal.indicator(stlInterval(3,4),true,false));

res = true;

% ------------------------------ END OF CODE ------------------------------
