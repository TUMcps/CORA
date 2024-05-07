function res = test_cora_warnings_enabled
% test_cora_warnings_enabled - tests the macro CORA_WARNINGS_ENABLED
%
% Syntax:
%    res = test_cora_warnings_enabled()
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
% See also: DEPRECATED_WARNINGS_ENABLED

% Authors:       Tobias Ladner
% Written:       14-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if all warnings are enabled
resvec = [
    % all warnings
    CORA_WARNINGS_ENABLED;
    CORA_WARNINGS_ENABLED('CORA:all');
    % folder specific warnings
    CORA_WARNINGS_ENABLED('CORA:app');
    CORA_WARNINGS_ENABLED('CORA:contDynamics');
    CORA_WARNINGS_ENABLED('CORA:contSet');
    CORA_WARNINGS_ENABLED('CORA:converter');
    CORA_WARNINGS_ENABLED('CORA:discDynamics');
    CORA_WARNINGS_ENABLED('CORA:examples');
    CORA_WARNINGS_ENABLED('CORA:global');
    CORA_WARNINGS_ENABLED('CORA:hybridDynamics');
    CORA_WARNINGS_ENABLED('CORA:manual');
    CORA_WARNINGS_ENABLED('CORA:matrixSets');
    CORA_WARNINGS_ENABLED('CORA:models');
    CORA_WARNINGS_ENABLED('CORA:specification');
    % special warnings
    CORA_WARNINGS_ENABLED('CORA:nn');
    CORA_WARNINGS_ENABLED('CORA:plot');
    CORA_WARNINGS_ENABLED('CORA:solver');
    CORA_WARNINGS_ENABLED('CORA:redundant');
    CORA_WARNINGS_ENABLED('CORA:deprecated');
]';

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
