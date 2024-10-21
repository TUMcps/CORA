function res = test_kleeneSignal_toBoolSignal
% test_kleeneSignal_toBoolSignal - unit test function of toBoolSignal
%
% Syntax:
%    res = test_kleeneSignal_toBoolSignal
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
boolSig = pointSegmentSignal([0,1,2,2.5],[true,true,false,true,true,false,false,true]);
unkInterval = stlInterval(1.5,2.25,true,false);
isUnknown = pointSegmentSignal.indicator(unkInterval,true,false);
sig1 = kleeneSignal(boolSig,isUnknown);

% no conversion of unknowns
try
    sig1.toBoolSignal();
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:notDefined')
        rethrow(ME)
    end
end

% conversion of unknowns to true
expected = boolSig.set(unkInterval,true);
assert(sig1.toBoolSignal(true) == expected);

% conversion of unknowns to false
expected = boolSig.set(unkInterval,false);
assert(sig1.toBoolSignal(false) == expected);

% signal without unknowns is never affected by the conversion
sig2 = kleeneSignal(boolSig,pointSegmentSignal.indicator(stlInterval(),false,false));
assert(sig2.toBoolSignal() == boolSig);
assert(sig2.toBoolSignal(true) == boolSig);
assert(sig2.toBoolSignal(false) == boolSig);

res = true;

% ------------------------------ END OF CODE ------------------------------
