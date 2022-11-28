function completed = example_converter_powerSystem2cora_all()
% example_converter_powerSystem2cora_all - example for creating all 
%    specified power system benchmarks; this function requires MATPOWER.
%
% Syntax:  
%    completed = example_converter_powerSystem2cora_all()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%    [1] M. Althoff, "Formal and Compositional Analysis of Power Systems 
%        using Reachable Sets", IEEE Transactions on Power Systems 29 (5), 
%        2014, 2270-2280

% Author:       Matthias Althoff
% Written:      14-April-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% list of MATPOWER models
models_MATPOWER = {'case4gs', 'case5', 'case6ww', 'case9', 'case14', ...
    'case24_ieee_rts', 'case30', 'case39', 'case57', 'case118', ...
    'case145', ...
    'case10ba', 'case12da', 'case15da', 'case15nbr', ...
    'case16am', 'case16ci', 'case17me', 'case18nbr', ...
    'case22', 'case28da', 'case33bw', 'case33mg', 'case34sa', ...
    'case38si', 'case51ga', 'case51he', 'case69', 'case70da', ...
    'case74ds', 'case85', 'case94pi', 'case118zh', 'case136ma', ...
    'case141'};


% set path
path = [CORAROOT filesep 'converter' filesep 'powerSystem2cora' filesep 'cases'];

%% create CORA models
for i = 1:length(models_MATPOWER)
    % get model
    model = models_MATPOWER{i}
    % load MATPOWER model
    loadPowerSystemCase(model,['MATPOWER_',model]);
    % load created power system description
    load([path filesep ['MATPOWER_',model]]);
    % convert to CORA
    powerSystem2cora(['MATPOWER_',model])
end

% example completed
completed = true;

%------------- END OF CODE --------------