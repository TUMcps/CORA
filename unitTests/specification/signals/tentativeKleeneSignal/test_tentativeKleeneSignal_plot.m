function res = test_tentativeKleeneSignal_plot
% test_tentativeKleeneSignal_plot - unit test function of plot
%
% Syntax:
%    res = test_tentativeKleeneSignal_plot
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
canBeTrueSig = pointSegmentSignal.indicator(stlInterval(0,2),true,false);
canBeFalseSig = pointSegmentSignal.indicator(stlInterval(1,3),true,false);
sig = tentativeKleeneSignal(canBeTrueSig,canBeFalseSig);

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
