function res = test_spectraShadow_printSet
% test_spectraShadow_printSet - unit test function of printSet
%
% Syntax:
%    res = test_spectraShadow_printSet
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

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test empty
SpS = spectraShadow.empty(2);

printSet(SpS)
printSet(SpS,'high')
printSet(SpS,'high',true)
printSet(SpS,'high',false)

% test normal set
A0 = eye(3);
A1 = [0 1 0;1 0 0;0 0 0];
A2 = [0 0 1;0 0 0;1 0 0];
SpS = spectraShadow([A0,A1,A2],[-1;2],3*eye(2));

printSet(SpS)
printSet(SpS,'high')
printSet(SpS,'high',true)
printSet(SpS,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
