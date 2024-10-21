function options = validateRLoptions(options)
% validateRLoptions - checks input and sets default values for the
%    options.rl struct for the agentRL function.
%
% Syntax:
%    options = nnHelper.validateRLoptions(options)
%
% Inputs:
%    options - struct (see agentRL)
%
% Outputs:
%    options - updated options
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: agentRL

% Authors:       Lukas Koller
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
persistent defaultRLfields
if isempty(defaultRLfields)
    defaultRLfields = {
        'gamma', .99;        
        'tau', .005;
        'expNoise', .2;
        'expNoiseTarget', .2;
        'expNoiseType', 'OU';
        'expDecayFactor', 1;
        'batchsize', 64;
        'buffersize', 1e6;
        'noise', .1;
        'earlyStop', inf;
        'printFreq', 50;
        'visRate', 50;
    };
end

persistent defaultActorTrainFields
if isempty(defaultActorTrainFields)
    defaultActorTrainFields = {
        'optim', nnAdamOptimizer(1e-4,.9,.999,1e-8,0);
        'eta', 0.1;
        'omega', 0.5;
        'exact_backprop', false;
        'zonotope_weight_update', 'outer_product';
    };
end

persistent defaultCriticTrainFields
if isempty(defaultCriticTrainFields)
    defaultCriticTrainFields = {
        'optim', nnAdamOptimizer(1e-3,.9,.999,1e-8,1e-2);
        'eta', 0.01;
        'exact_backprop', false;
        'zonotope_weight_update', 'outer_product';
    };
end

persistent defaultAdvOps
if isempty(defaultAdvOps)
    defaultAdvOps = {
       'numSamples', 200;
       'alpha', 4;
       'beta', 4;
    };
end

% Validate neural network options. Do not set training options.
options = nnHelper.validateNNoptions(options,false);

% Set default RL options.
if ~isfield(options,'rl')
    options.rl = struct;
end
options.rl = nnHelper.setDefaultFields(options.rl,defaultRLfields);

% Set default actor training options.
if ~isfield(options.rl,'actor') || ~isfield(options.rl.actor,'nn') ...
        || ~isfield(options.rl.actor.nn,'train')
    options.rl.actor.nn.train = struct;
end
options.rl.actor.nn.train = ...
    nnHelper.setDefaultFields(options.rl.actor.nn.train,defaultActorTrainFields);
% Validate actor options. Set training options.
options.rl.actor = ...
    nnHelper.validateNNoptions(options.rl.actor,true);

% Set default critic training options.
if ~isfield(options.rl,'critic') || ~isfield(options.rl.critic,'nn') ...
        || ~isfield(options.rl.critic.nn,'train')
    options.rl.critic.nn.train = struct;
end
options.rl.critic.nn.train = ...
    nnHelper.setDefaultFields(options.rl.critic.nn.train,defaultCriticTrainFields);
% Validate critic options. Set training options.
options.rl.critic = ...
    nnHelper.validateNNoptions(options.rl.critic,true);

% Validate adversarial attack options.
if ~isfield(options.rl.actor.nn.train,'advOps')
    options.rl.actor.nn.train.advOps = struct;
end
options.rl.actor.nn.train.advOps = ...
    nnHelper.setDefaultFields(options.rl.actor.nn.train.advOps,defaultAdvOps);

% Check rl fields
if CHECKS_ENABLED
    structName = inputname(1);
    % Check rl fields
    aux_checkFieldNumericDefInterval(options.rl,'gamma',interval(0,1),structName)
    aux_checkFieldNumericDefInterval(options.rl,'tau',interval(0,1),structName)
    aux_checkFieldNumericDefInterval(options.rl,'expNoise',interval(0,inf),structName)
    aux_checkFieldNumericDefInterval(options.rl,'expNoiseTarget',interval(0,inf),structName)
    aux_checkFieldStr(options.rl, 'expNoiseType', {'OU','gaussian'}, structName);
    aux_checkFieldNumericDefInterval(options.rl,'expDecayFactor',interval(-1,1),structName)
    aux_checkFieldNumericDefInterval(options.rl,'batchsize',interval(0,inf),structName)
    aux_checkFieldNumericDefInterval(options.rl,'noise',interval(0,inf),structName)

    % Check actor fields
    aux_checkFieldNumericDefInterval(options.rl.actor.nn.train,'eta',interval(0,inf),structName)
    aux_checkFieldNumericDefInterval(options.rl.actor.nn.train,'omega',interval(0,1),structName)

    % Check critic fields
    aux_checkFieldNumericDefInterval(options.rl.critic.nn.train,'eta',interval(0,inf),structName)

    % Validate the train options 
    if strcmp(options.rl.critic.nn.train.method,'set')
        if ~strcmp(options.rl.actor.nn.train.method,'set')
            throw(CORAerror('CORA:wrongFieldValue', ...
                'critic.nn.train.method', "set"))
        end
    end
    aux_checkFieldStr(options.rl.critic.nn.train,'method',{'point','set'},structName);
end

end


% Auxiliary functions -----------------------------------------------------

function aux_checkFieldStr(optionsnn, field, admissibleValues, structName)
    fieldValue = optionsnn.(field);
    if ~(isa(fieldValue, 'string') || isa(fieldValue, 'char')) || ...
        ~ismember(fieldValue, admissibleValues)
        throw(CORAerror('CORA:wrongFieldValue', ...
            aux_getName(structName, field), admissibleValues))
    end
end

function aux_checkFieldNumericDefInterval(optionsrl, field, I, structName)
    if ~contains_(I,optionsrl.(field),'exact',eps)
        throw(CORAerror('CORA:outOfDomain', ...
            aux_getName(structName, field), "ValidDomain", I))
    end
end

function msg = aux_getName(structName, field)
    msg = sprintf("%s.nn.%s", structName, field);
end

% ------------------------------ END OF CODE ------------------------------
