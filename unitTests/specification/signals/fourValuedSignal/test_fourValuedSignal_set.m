function res = test_fourValuedSignal_set
% test_fourValuedSignal_set - unit test function of set
%
% Syntax:
%    res = test_fourValuedSignal_set
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

res = true;

tt = fourValued.True;
uu = fourValued.Unknown;
ff = fourValued.False;
ii = fourValued.Inconclusive;

% signals
boolSig = pointSegmentSignal([0,1,2,2.5],[true,true,false,true,true,false,false,true]);
isUnknown = pointSegmentSignal.indicator(stlInterval(1.5,2.25,true,false),true,false);
isInconclusive = pointSegmentSignal.indicator(stlInterval(1.75,2.1,false,true),true,false);
kleeneSig = kleeneSignal(boolSig,isUnknown);
sig = fourValuedSignal(kleeneSig,isInconclusive);

% set to inconclusive
actual = sig.set(stlInterval(0,0.5),ii);
assert(actual.isInconclusive == isInconclusive.set(stlInterval(0,0.5),true));
actual = sig.set(stlInterval(0.5,1.75,false),ii);
assert(actual.isInconclusive == isInconclusive.set(stlInterval(0.5,1.75,false),true));

% set to true
actual = sig.set(stlInterval(2.3,2.5),tt);
assert(actual.isInconclusive == isInconclusive);
assert(actual.signal == kleeneSig.set(stlInterval(2.3,2.5),kleene.True));
actual = sig.set(stlInterval(1.75,2),tt);
assert(actual.isInconclusive == isInconclusive.set(stlInterval(1.75,2),false));
assert(actual.signal == kleeneSig.set(stlInterval(1.75,2),kleene.True));

% set to unknown
actual = sig.set(stlInterval(0,1.5),uu);
assert(actual.isInconclusive == isInconclusive);
assert(actual.signal == kleeneSig.set(stlInterval(0,1.5),kleene.Unknown));
actual = sig.set(stlInterval(1.6,2.2),uu);
assert(actual.isInconclusive == isInconclusive.set(stlInterval(1.6,2.2),false));
assert(actual.signal == kleeneSig.set(stlInterval(1.6,2.2),kleene.Unknown));

% set to false
actual = sig.set(stlInterval(0,1.25),ff);
assert(actual.isInconclusive == isInconclusive);
assert(actual.signal == kleeneSig.set(stlInterval(0,1.25),kleene.False));
actual = sig.set(stlInterval(1.25,1.75),ff);
assert(actual.isInconclusive == isInconclusive.set(stlInterval(1.25,1.75),false));
assert(actual.signal == kleeneSig.set(stlInterval(1.25,1.75),kleene.False));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
