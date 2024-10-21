function res = test_fourValuedSignal_toKleeneSignal
% test_fourValuedSignal_toKleeneSignal - unit test function of toKleeneSignal
%
% Syntax:
%    res = test_fourValuedSignal_toKleeneSignal
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
boolSig = pointSegmentSignal([0,1,2,2.5],[true,true,false,true,true,false,false,true]);
isUnknown = pointSegmentSignal.indicator(stlInterval(1.5,2.25,true,false),true,false);
kleeneSig = kleeneSignal(boolSig,isUnknown);
incInterval = stlInterval(1.75,2.1,false,true);
isInconclusive = pointSegmentSignal.indicator(incInterval,true,false);
sig1 = fourValuedSignal(kleeneSig,isInconclusive);

% no conversion of inconclusives
try
    sig1.toKleeneSignal();
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:notDefined')
        rethrow(ME)
    end
end

% conversion of inconclusives to true
expected = kleeneSig.set(incInterval,kleene.True);
assert(sig1.toKleeneSignal(kleene.True) == expected);

% conversion of inconclusives to unknown
expected = kleeneSig.set(incInterval,kleene.Unknown);
assert(sig1.toKleeneSignal(kleene.Unknown) == expected);

% conversion of inconclusives to false
expected = kleeneSig.set(incInterval,kleene.False);
assert(sig1.toKleeneSignal(kleene.False) == expected);

% signal without inconclusives is never affected by the conversion
sig2 = fourValuedSignal(kleeneSig,pointSegmentSignal.indicator(stlInterval(),false,false));
assert(sig2.toKleeneSignal() == kleeneSig);
assert(sig2.toKleeneSignal(kleene.True) == kleeneSig);
assert(sig2.toKleeneSignal(kleene.Unknown) == kleeneSig);
assert(sig2.toKleeneSignal(kleene.False) == kleeneSig);

res = true;

% ------------------------------ END OF CODE ------------------------------
