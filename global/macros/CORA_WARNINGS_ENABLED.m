function res = CORA_WARNINGS_ENABLED(identifier)
% CORA_WARNINGS_ENABLED - specifies if a CORA warning should be shown
%
% Syntax:
%    res = CORA_WARNINGS_ENABLED()
%
% Inputs:
%    identifier - char
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORAwarning

% Authors:       Tobias Ladner
% Written:       14-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 1
    identifier = 'CORA:all';
end

% check if enabled
switch identifier

    case 'CORA:all'
        % controls all CORA warnings
        res = true;

        % general folder warnings -----------------------------------------

    case 'CORA:app'
        % controls general CORA warnings in ./app
        res = true;

    case 'CORA:contDynamics'
        % controls general CORA warnings in ./contDynamics
        res = true;

    case 'CORA:contSet'
        % controls general CORA warnings in ./contSet
        res = true;

    case 'CORA:converter'
        % controls general CORA warnings in ./converter
        res = true;

    case 'CORA:discDynamics'
        % controls general CORA warnings in ./discDynamics
        res = true;

    case 'CORA:examples'
        % controls general CORA warnings in ./examples
        res = true;

    case 'CORA:global'
        % controls general CORA warnings in ./global
        res = true;

    case 'CORA:hybridDynamics'
        % controls general CORA warnings in ./hybridDynamics
        res = true;

    case 'CORA:manual'
        % controls general CORA warnings in ./manual
        res = true;

    case 'CORA:matrixSets'
        % controls general CORA warnings in ./matrixSets
        res = true;

    case 'CORA:models'
        % controls general CORA warnings in ./models
        res = true;

    case 'CORA:specification'
        % controls general CORA warnings in ./specification
        res = true;

        % special warnings ------------------------------------------------

    case 'CORA:nn'
        % controls CORA warnings regarding neural network verification
        res = true;

    case 'CORA:solver'
        % controls CORA warnings regarding solvers
        res = true;

    case 'CORA:plot'
        % controls CORA warnings regarding plotting
        res = true;
    
    case 'CORA:deprecated'
        % controls deprecation warnings
        res = true;
    
    case 'CORA:redundant'
        % controls redundancy warnings
        res = true;

    otherwise
        throw(CORAerror('CORA:wrongValue','first',sprintf('Unknown identifier ''%s''',identifier)))
end


% ------------------------------ END OF CODE ------------------------------
