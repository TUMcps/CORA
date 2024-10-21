function res = test_kleeneSignal_set
% test_kleeneSignal_set - unit test function of set
%
% Syntax:
%    res = test_kleeneSignal_set
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

res = true;

tt = kleene.True;
uu = kleene.Unknown;
ff = kleene.False;

% signals
boolSig = pointSegmentSignal([0,1,2,2.5],[true,true,false,true,true,false,false,true]);
isUnknown = pointSegmentSignal.indicator(stlInterval(1.5,2.25,true,false),true,false);
sig = kleeneSignal(boolSig,isUnknown);

% set to unknown
actual = sig.set(stlInterval(0,0.5),uu);
assert(actual.isUnknown == (isUnknown | pointSegmentSignal.indicator(stlInterval(0,0.5),true,false)));
actual = sig.set(stlInterval(0.5,1.75,false),uu);
assert(actual.isUnknown == (isUnknown | pointSegmentSignal.indicator(stlInterval(0.5,1.75,false),true,false)));

% set to true
actual = sig.set(stlInterval(2.3,2.5),tt);
assert(actual.isUnknown == isUnknown);
assert(actual.signal == (boolSig | pointSegmentSignal.indicator(stlInterval(2.3,2.5),true,false)));
actual = sig.set(stlInterval(1.75,2),tt);
assert(actual.isUnknown == (isUnknown & pointSegmentSignal.indicator(stlInterval(1.75,2),false,true)));
assert(actual.signal == (boolSig | pointSegmentSignal.indicator(stlInterval(1.75,2),true,false)));

% set to false
actual = sig.set(stlInterval(0,1.25),ff);
assert(actual.isUnknown == isUnknown);
assert(actual.signal == (boolSig & pointSegmentSignal.indicator(stlInterval(0,1.25),false,true)));
actual = sig.set(stlInterval(1.25,1.75),ff);
assert(actual.isUnknown == (isUnknown & pointSegmentSignal.indicator(stlInterval(1.25,1.75),false,true)));
assert(actual.signal == (boolSig & pointSegmentSignal.indicator(stlInterval(1.25,1.75),false,true)));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
