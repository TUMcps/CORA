function res = test_pointSegmentSignal_plot
% test_pointSegmentSignal_plot - unit test function of plot
%
% Syntax:
%    res = test_pointSegmentSignal_plot
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
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tt = true;
ff = false;
sig = pointSegmentSignal([0,1,2,2.5], [tt,tt,ff,tt,tt,ff,ff,tt]);

try
    figure;
    sig.plot();
    han = sig.plot(); %#ok<NASGU>
    close;
catch
    close;
    rethrow(ME)
end

res = true;

% ------------------------------ END OF CODE ------------------------------
