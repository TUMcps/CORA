function pass = test_spaceex2cora_parseConfig()
% test_spaceex2cora_parseConfig - example for parsing of SpaceEx
%                                 configuration file
%
% Syntax:
%    test_spaceex2cora_parseConfig
%
% Inputs:
%    -
%
% Outputs:
%    pass - boolean

% Authors:       Maximilian Perschl
% Written:       23-September-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pass = true;

% Since parsing of initial and forbidden conditions (which are the main
% tasks of the configuration-converter) are check within their own seperate
% tests, this function checks if the file reading and data transfer of the
% function "parseSpaceExConfig" is working as intended.

% Carry out conversion, with added configuration file
spaceex2cora('hybrid_bball.xml',[],[],[],[],'bball.cfg');
% Load conversion results
load('hybrid_bball_config.mat','configParams','configSpecs','spec_mapping');
% Check conversion results for correctness
% parameters
expectedParams.Tfinal = 40;
expectedParams.Tsample = 0.1;
expectedParams.R0 = interval([10;0],[10.2;0]);
expectedParams.u = interval.empty(1);
if ~isequal(expectedParams,configParams)
    pass = false; return
end

% specifications
expectedSpecs = specification(halfspace.empty(2));
A = [-1 0]; b = -100;
expectedSpecs = add(expectedSpecs,specification(polytope(A,b)));
% Specifciations are compared for number of specifications, type, and
% corresponding set
if length(expectedSpecs) ~= length(configSpecs)
    pass = false; return
else
    % Start at 2 since we ignore the first specification which has an empty
    % set
    for i = 2:length(expectedSpecs)
        if ~(contains(expectedSpecs(i).set,configSpecs(i).set) && ...
                contains(configSpecs(i).set,expectedSpecs(i).set) && ...
                isequal(configSpecs(i).type,expectedSpecs(i).type))
            pass = false; return
        end
    end
end
% specification mapping
if ~isequal(spec_mapping,{})
    pass = false;
end

end

% ------------------------------ END OF CODE ------------------------------
