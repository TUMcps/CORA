function res = test_tentativeKleeneSignal_emptySignal
% test_tentativeKleeneSignal_emptySignal - unit test function of emptySignal
%
% Syntax:
%    res = test_tentativeKleeneSignal_emptySignal
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
 
% assume true
res = true;

% signals
empty = tentativeKleeneSignal.emptySignal();
blank = pointSegmentSignal.indicator(stlInterval(),false,false);

assert(empty.canBeTrue == blank && empty.canBeFalse == blank);

% ------------------------------ END OF CODE ------------------------------
