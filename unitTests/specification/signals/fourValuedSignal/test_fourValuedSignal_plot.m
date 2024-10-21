function res = test_fourValuedSignal_plot
% test_fourValuedSignal_plot - unit test function of plot
%
% Syntax:
%    res = test_fourValuedSignal_plot
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

boolSig = pointSegmentSignal([0,1,2,2.5],[true,true,false,true,true,false,false,true]);
isUnknown = pointSegmentSignal.indicator(stlInterval(1.5,2.25,true,false),true,false);
isInconclusive = pointSegmentSignal.indicator(stlInterval(1.75,2.1,false,true),true,false);
sig = fourValuedSignal(kleeneSignal(boolSig,isUnknown),isInconclusive);

try
    figure;
    sig.plot();
    han = sig.plot(); %#ok<NASGU>
    close;
catch ME
    close;
    rethrow(ME)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
