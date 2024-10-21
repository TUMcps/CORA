function res = test_kleeneSignal_plot
% test_kleeneSignal_plot - unit test function of plot
%
% Syntax:
%    res = test_kleeneSignal_plot
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

tt = true;
ff = false;
boolSig = pointSegmentSignal([0,1,2,2.5],[tt,tt,ff,tt,tt,ff,ff,tt]);
isUnknown = pointSegmentSignal.indicator(stlInterval(1.5,2.25,tt,ff),tt,ff);
sig = kleeneSignal(boolSig,isUnknown);

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
