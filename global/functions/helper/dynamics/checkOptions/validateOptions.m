function options = validateOptions(sys,func,params,options,varargin)
% validateOptions - validates model parameters and algorithm parameters
%    for given function
%
% Syntax:
%    options = validateOptions(sys,func,params,options)
%    options = validateOptions(sys,func,params,options,internalCall)
%
% Inputs:
%    sys - object of contDynamics class
%    func - function for which parameters should be validated
%    params - model parameters (e.g., initial set)
%    options - algorithm parameters (e.g., time step size)
%    internalCall - (optional) true/false whether validation has happened
%                   before at a higher system level; allows to skip
%                   validation and only do post-processing
%
% Outputs:
%    options -  merged parameters and options for further computations
%
% Other m-files required:
%    initErrorCodex.m - initializes error identifiers used to describe
%                       errors that occur during the validation process
%    getDefaultValue.m - stores the default values for all parameters that
%                        have default values but are not set by the user
%    getErrorMessage.m - returns the error message corresponding to the
%                        identifier defined in initErrorCodex.m (only
%                        relevant if any parameter is set incorrectly)
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       25-January-2021
% Last update:   26-January-2021
%                05-October-2023 (TL, simplified config files)
% Last revision: 05-January-2023 (MW, extend description of code, no more errors but text on console)
%                19-June-2023 (MW, adapt to new configfile syntax)

% ------------------------------ BEGIN CODE -------------------------------

% set default value for internal call
internalCall = setDefaultValues({false},varargin);

% get file name of configuration file
if ~internalCall
    [fileFound,configfile,cfgIdx] = aux_getConfigfile(sys,func);
end

if ~internalCall && ~fileFound
    
    % if no configuration file exists, we notify the user that this is the
    % case (occurs most likely during development) and that parameters and
    % options are not being checked

    disp("No configuration file to execute model/algorithm parameter validation.");
    disp("-> model/algorithm parameters are not checked!");
    
else
    % params/options validation in multiple steps:

    % we first check the model parameters, then the algorithm parameters;
    % this is because the model parameters do not depend on the setting of
    % the algorithm parameters (cases such as the length of an input signal
    % and the time step size exist, but then it is deemed to be the fault
    % of the algorithm parameter (time step size) which does not match the
    % length of the input signal
    % the check of either model or algorithm parameters does not include
    % any rewriting of the provided values (this is done in step 4.)
    
    if ~internalCall
    % 2. validation of model parameters (params)
    params = aux_checkFullList(sys,func,params,options,configfile,'params');
    
    % 3. validation of algorithm parameters (options)
    options = aux_checkFullList(sys,func,params,options,configfile,'options');
    end
    
    % 4. post-processing: internal rewriting of the checked model/algorithm
    % parameters to allow for a simpler syntax in the reach/simulate/...
    % functions
    [params, options] = postProcessing(sys,func,params,options);
    
    % 5. since the split in params and options historically occurred after
    % many of the functions had already been written, we unify the structs
    % for params and options and work with only options in the code
    options = params2options(params,options);
end

end


% Auxiliary functions -----------------------------------------------------

function [fileFound,configfile,cfgIdx] = aux_getConfigfile(sys,func)
% sys - object of contDynamics class
% func - function for which parameters/options should be validated
%
% fileFound - true/false whether a configuration file was found
%             (false has special handling for development)
% configfile - list of names of configuration files
% cfgIdx - list of indices: for all contDynamics systems just 0, otherwise
%          0 for root component of hybridAutomaton/parallelHybridAutomaton
%          and indices of locations with additional configuration file
%          note: for parallel hybrid automata, we use linear indices, e.g.,
%          '5' for location 3 of subcomponent 2, when subcomponent 1 has 2
%          locations, etc.
%
% in this auxiliary function, we check whether a configuration file exists
% for the given pair of system and function; in the case of hybrid
% dynamics, we also look through all locations and return the corresponding
% configuration files (e.g., for linearSys/reach for a location with linear
% dynamics in a hybrid automaton)

currfolder = mfilename('fullpath');
currfolder = currfolder(1:end-length(mfilename));

% read name of contDynamics class
classname = class(sys);

% rename classname to 'contDynamics' as there is only one configuration
% file for all classes inheriting from contDynamics
if isa(sys,'contDynamics') && any(contains(func,{'simulate','observe'}))
    classname = 'contDynamics';
end

% obtain name of configfile for given class and function
configfile = {['config_' classname '_' func]};
cfgIdx = 0;

% check if configfile exists
fileFound = isfile([currfolder filesep 'configfiles' filesep configfile{1} '.m']);


% if the system is hybrid, include configuration files of locations (of
% subcomponents), too: to this end, we go through all locations (of all
% subcomponents) and evaluate the corresponding configuration file to
% obtain the list of parameters and options
if (isa(sys,'hybridAutomaton') || isa(sys,'parallelHybridAutomaton')) ...
        && ~strcmp(func,'simulate')
    % simulateRandom -> simulateStandard (corresponding func in contDynamics)

    if isa(sys,'hybridAutomaton')
        % init list of additional configuration files
        addConfigfiles = cell(length(sys.location),1);

        % loop over all locations
        for i=1:length(sys.location)
            % read out contDynamics class
            classname = class(sys.location(i).contDynamics);
            if contains(func,'simulate')
                % for simulateRandom|simulate just contDynamics
                classname = 'contDynamics';
            end
            % append name of configuration file
            addConfigfiles{i} = ['config_' classname '_' func];
        end

    elseif isa(sys,'parallelHybridAutomaton')

        % total number of individual locations
        numLocs = arrayfun(@(x) length(x.location),sys.components,'UniformOutput',true);
        % init list of configuration files
        addConfigfiles = cell(sum(numLocs),1);

        % loop over all components
        for i=1:length(sys.components)

            % loop over all locations in i-th component
            for j=1:length(sys.components(i).location)

                % read out contDynamics class of i/j-th location
                classname = class(sys.components(i).location(j).contDynamics);
                if contains(func,'simulate')
                    % for simulateNormal|simulate just contDynamics
                    classname = 'contDynamics';
                end
                addConfigfiles{sum(numLocs(1:i-1))+j,1} = ...
                    ['config_' classname '_' func];
            end
        end
    end

    % remove all duplicates, keep indices of remaining names
    [addConfigfiles,addCfgIdx] = unique(addConfigfiles);

    % append to list of configuration files
    configfile = [configfile; addConfigfiles];
    cfgIdx = [cfgIdx; addCfgIdx];

end

end

function blueprint = aux_readConfigfiles(sys,params,options,configfile,listname)
% sys           ... object of some contDynamics class
% params        ... list of model parameters
% options       ... list of algorithmic parameters
% configfile    ... name of configuration file (for hybrid systems: list of
%                   unique names)
% cfgIdx        ... list of indices for locations with unique configuration
%                   file (0 for purely continuous systems)
% listname      ... 'params', 'options', for error messages

% loop over all configuration files
for i=1:length(configfile)
    if strcmp(listname,'params')
        [addBlueprint,~] = eval([configfile{i} '(sys,params,options)']);
    elseif strcmp(listname,'options')
        [~,addBlueprint] = eval([configfile{i} '(sys,params,options)']);
    end
    if i == 1
        blueprint = addBlueprint;
    else
        idxNew = ~ismember({addBlueprint.name}',{blueprint.name}');
        blueprint = [blueprint; addBlueprint(idxNew)];
    end
end

% read condition function
for i=1:length(blueprint)
    blueprint(i).condfun = getCondfunDynParameter(blueprint(i).name,listname);
end

end

function list = aux_checkFullList(sys,func,params,options,configfile,listname)
% sys           ... object of some contDynamics class
% func          ... name of checked function, i.e. 'reach'
% params        ... list of model parameters
% options       ... list of algorithmic parameters
% configfile    ... name of configuration file (for hybrid systems: list of
%                   unique names)
% cfgIdx        ... list of indices for locations with unique configuration
%                   file (0 for purely continuous systems)
% listname      ... 'params', 'options', for error messages

% in this function, we check the provided list of model parameters (params)
% or algorithm parameters (options); due to the occurrence of conditional
% parameters (i.e., parameters that are mandatory only if another parameter
% has a certain value), the sequence of steps is somewhat cumbersome:
% 1. set non-conditional defaults
%    for all parameters,
%       for which there exists a default value
%       and which do not depend on another parameter,
%    we check whether the user has provided some value for this parameter;
%    if not, then the default value is assigned to this parameter
%
% 2. check for redundant and conditional parameters
%    we now check which of the user-provided model/algorithm parameters are
%       redundant (not expected to exist)
%       and conditional (depending on other parameters)
%    to exclude them for the subsequent first round of parameter checks
%
% 3. after all non-conditional parameters have been checked, we can now
%    know whether the conditions for the conditional parameters are
%    fulfilled and checked those parameters as well
%
% 4. all redundancies are printed (if VALIDATEOPTIONS_ERRORS = false)

% read the configuration file(s)
blueprint = aux_readConfigfiles(sys,params,options,configfile,listname);


% 1. set non-conditional defaults -----------------------------------------

list = aux_setMissingDefaults(params,options,blueprint,sys,listname);

% update configuration file
if strcmp(listname,'params')
    blueprint = aux_readConfigfiles(sys,list,options,configfile,listname);
elseif strcmp(listname,'options')
    blueprint = aux_readConfigfiles(sys,params,list,configfile,listname);
end


% 2. check irredundant and non-conditional parameters ---------------------

% read all fields (special read out when field is a struct)
allfields = aux_readAllFields(list);

% redundant indices
res.redIdxInUserList = false(length(allfields),1);

% conditional indices
condIdx = false(length(allfields),1);

% go through all fields, see which are conditional or redundant
for i=1:length(allfields)

    % find idx of current field in blueprint
    idx = find(ismember({blueprint.name}',allfields{i}),1,'first');

    if isempty(idx)
        % field not in blueprint -> param / option is redundant
        res.redIdxInUserList(i) = true;
    elseif ~isempty(blueprint(idx).condfun)
        % field not in blueprint -> param / option is redundant
        condIdx(i) = true;
    end
end

% skip conditional and redundant indices
skipIdx = condIdx | res.redIdxInUserList;

% check non-conditional and non-redundant params / options
[res.mandMissing,res.failedChecks] = aux_validateList(list,blueprint,sys,func,params,options,listname,skipIdx);


% 3. check remaining parameters -------------------------------------------

% check conditions for conditional params / options
% ... remove 'cond-' from status
if strcmp(listname,'params')
    [blueprint,condBPIdx,condIdx] = aux_checkCondition(blueprint,func,listname,sys,list,options);
elseif strcmp(listname,'options')
    [blueprint,condBPIdx,condIdx] = aux_checkCondition(blueprint,func,listname,sys,params,list);
end

if any(condBPIdx)
    % set conditional defaults (now also listed in blueprint as 'default')
    if strcmp(listname,'params')
        list = aux_setMissingDefaults(list,options,blueprint,sys,listname);
    elseif strcmp(listname,'options')
        list = aux_setMissingDefaults(params,list,blueprint,sys,listname);
    end
    
    % another update of blueprint would be required if any of the
    % conditional defaults were set automatically to their default value
    % and they would be in a check function of another conditional param / option  
    
    % unify condIdx
    condIdx = [condIdx; true(length(aux_readAllFields(list))-length(condIdx),1)];
    
    % remaining: check conditional params / options (also cond-default...)
    [mandMissing_,failedChecks_] = aux_validateList(list,blueprint,sys,func,params,options,listname,~condIdx);

    % append to previous result struct
    res.mandMissing = res.mandMissing | mandMissing_;
    res.failedChecks = [res.failedChecks; failedChecks_];
end


% 4. print information about input argument validation --------------------
% depends on value of macro VALIDATEOPTIONS_ERRORS:
% - true:  warning for redundant params/options, since they are deleted from
%          the respective struct
% - false: parameters which are not mentioned in the configfile, missing
%          mandatory parameters, and failed checks

if strcmp(listname,'params')
    aux_printValidationResult(list,options,blueprint,sys,func,listname,res);
elseif strcmp(listname,'options')
    aux_printValidationResult(params,list,blueprint,sys,func,listname,res);
end

end


% Auxiliary functions -----------------------------------------------------

function list = aux_setMissingDefaults(params,options,blueprint,sys,listname)
% params    ... params struct
% options   ... options struct
% blueprint ... paramsList / optionsList (from config file)
% sys       ... object of some hybrid or contDynamics class
% listname  ... 'params' or 'options'

% here, we check whether model/algorithm parameters with default values are
% assigned a value by the user:
%    if yes, proceed to next parameter
%    if not, assign the default value
% note: the default values are defined in the function getDefaultValue.m

% select struct
if strcmp(listname,'params')
    list = params;
elseif strcmp(listname,'options')
    list = options;
end

% indices in the blueprint of parameters with default values
% (here ,we skip all parameters which are conditional, i.e., whose status
% depends on the value of another parameter / which has a condfun)
idxDefault = ismember({blueprint.status}','default') & ...
    cellfun(@(x)isempty(x),{blueprint.condfun}');

% default fields
defFields = {blueprint(idxDefault).name}';

% names of user-defined fields (including special handling for parameters
% which are themselves structs)
fieldnames = aux_readAllFields(list);

% loop over default params / options
for i=1:length(defFields)

    % assign default value by reading from list of default values
    if ~any(ismember(fieldnames,defFields{i}))

        % assign default value by reading from list of default values
        if ~contains(defFields{i},'.')
            % standard (no struct) parameter

            % assign default value (full list in getDefaultValue.m)
            if strcmp(listname,'params')
                list.(defFields{i}) = ...
                    getDefaultValue(defFields{i},sys,list,options,listname);
            elseif strcmp(listname,'options')
                list.(defFields{i}) = ...
                    getDefaultValue(defFields{i},sys,params,list,listname);
            end
        else
            % parameter is a struct -> read out name of struct and field
            dotIdx = strfind(defFields{i},'.');
            firstname = defFields{i}(1:dotIdx-1);
            secondname = defFields{i}(dotIdx+1:end);

            % assign default value (full list in getDefaultValue.m)
            if strcmp(listname,'params')
                list.(firstname).(secondname) = ...
                    getDefaultValue(defFields{i},sys,list,options,listname);
            elseif strcmp(listname,'options')
                list.(firstname).(secondname) = ...
                    getDefaultValue(defFields{i},sys,params,list,listname);
            end
            
        end
    end

end

end

function [mandMissing,failedChecks] = aux_validateList(list,blueprint,sys,func,params,options,listname,skipIdx)
% list      ... params struct / options struct
% blueprint ... paramsList / optionsList (from config file)
% sys       ... object of some hybridDynamics or contDynamics class
% listname  ... {'params','options'} (for error messages)
% skipIdx   ... indices in list which are to be skipped
%               (redundant and conditional params / options)

% for each parameter in 'list', we check whether it is provided by the user
% (only if status = 'mandatory') and evaluate the check functions defined
% in the corresponding entry in 'blueprint' (from the configfile)
% -> wrongdoing is saved in output arguments

% read out fields from list (including special handling for structs)
allfields = aux_readAllFields(list);

% monitor which parameters defined as mandatory in 'blueprint' are missing
% in the user-provided set of parameters 'list'
mandMissing = false(length(blueprint),1);


% loop over all entries in 'blueprint'
for i=1:length(blueprint)

    % check status (and existence of condfun)
    if strcmp(blueprint(i).status,'mandatory') && isempty(blueprint(i).condfun)
        % i-th parameter in 'blueprint' is mandatory and does not depend on
        % the setting of another parameter (empty condfun)

        % check if the user has provided a value for this parameter
        if ~isfield(list,blueprint(i).name)

            if VALIDATEOPTIONS_ERRORS
                % mandatory parameter missing -> print error
                throw(CORAerror('CORA:specialError',...
                    sprintf('Error in %s check for %s object:\n  %s.%s is missing.',...
                    listname,class(sys),listname,blueprint(i).name)));
            else
                % save index of missing parameter
                mandMissing(i) = true;
            end
        end
    end
end

% init indices for which a check has failed (note: the '[]' is crucial for
% the concatenation later on!)
failedChecks = struct([]);

% go over all parameters (also the ones to which we have assigned the
% default value in aux_setMissingDefaults) in the order by which they
% are sorted in the configfile
for i=1:length(blueprint)
    
    % for i-th parameter in blueprint, find the corresponding index in
    % the user-provided list of parameters
    idx = find(ismember(allfields,blueprint(i).name),1,'first');
    
    % skip if already successful check (or to be skipped b/c cond)
    if ~isempty(idx) && ~skipIdx(idx) && ~mandMissing(i)
        % check field
        if strcmp(listname,'params')
            failedChecks_i = checkDynParameter(blueprint(i).name,sys,func,list,options,listname);
        elseif strcmp(listname,'options')
            failedChecks_i = checkDynParameter(blueprint(i).name,sys,func,params,list,listname);
        end
        failedChecks = [failedChecks,failedChecks_i];
    end
end

end

function [blueprint,condBPIdx,condIdx] = aux_checkCondition(blueprint,func,listname,sys,params,options)
% list      ... params struct / options struct
% blueprint ... paramsList / optionsList (from config file)

% this function returns the indices of the entries in 'list' and
% 'blueprint', where the stored parameters are conditional; additionally,
% we remove the condition function from the blueprint

% find out which indices in 'blueprint' are of parameters whose status
% depends on other parameters (have a non-empty condfun)
condBPIdx = false(length(blueprint),1);

% loop over all entries in 'blueprint'
for i=1:length(blueprint)

    if ~isempty(blueprint(i).condfun)
        condBPIdx(i) = true;

        condRes = blueprint(i).condfun(sys,func,params,options);

        % evaluate the condition function: if it is not fulfilled, the 
        % parameter's status is set to 'redundant'
        if ~condRes
            % set field 'status' to 'redundant'
            blueprint(i).status = 'redundant';
            condBPIdx(i) = false;
        end


    end

    if ~isempty(blueprint(i).condfun)
        % remove the condition function from the blueprint
        blueprint(i).condfun = [];
    end
end

% now, determine indices in 'list' which corresponding to the indices in
% 'blueprint' which are conditional parameters

% read out all fields from the user-provided list of parameters (including
% special handling for parameter which are themselves structs)
if strcmp(listname,'params')
    allfields = aux_readAllFields(params);
elseif strcmp(listname,'options')
    allfields = aux_readAllFields(options);
end

% init indices for conditional parameters (indexing 'list')
condIdx = false(length(allfields),1);

% loop over all entries in 'blueprint'
for i=1:length(condBPIdx)
    % check if i-th entry in 'blueprint' is a conditional parameter
    if condBPIdx(i)
        % corresponding index in 'list' for i-th entry in 'blueprint'
        idx = find(ismember(allfields,blueprint(i).name),1,'first');
        % set conditional index to true
        condIdx(idx) = true;
    end
end

end


function aux_printValidationResult(params,options,blueprint,sys,func,listname,res)
% params    ... params struct
% options   ... options struct
% blueprint ... paramsList / optionsList (from config file)
% sys       ... object of some hybrid or contDynamics class
% func      ... name of checked function, i.e. 'reach'
% listname  ... 'params' or 'options'
% res       ... struct containing the result of the parameter validation

% prints the result of the parameter validation in case infractions like
% 1.  missing mandatory parameters
% 2.  parameters fail check functions
% 3.  redundant parameters (not defined in configfile)
% are found (1 and 2 only if VALIDATEOPTIONS_ERRORS = false)

% no infractions found
if ~any(res.mandMissing) && isempty(res.failedChecks) ...
        && ~any(res.redIdxInUserList)
    return
end

headerPrinted = false;

% print missing mandatory values and failed checks
if ~VALIDATEOPTIONS_ERRORS

    % 1. missing mandatory values
    mandMissingTxt = ['''', strjoin({blueprint(res.mandMissing).name},', '), ''''];
    if ~strcmp(mandMissingTxt,"''")
        % print header
        headerPrinted = aux_printHeader(listname,headerPrinted);
        % print missing mandatory fields
        fprintf("  - missing mandatory fields: " + mandMissingTxt + "\n");
    end
    
    
    % 2. failed checks (with message)
    numFailedChecks = length(res.failedChecks);
    for i=1:numFailedChecks
        % print header (unless already printed)
        headerPrinted = aux_printHeader(listname,headerPrinted);

        % check how many functions failed for the same parameter
        fieldname = res.failedChecks(i).parameter;
        j = i;
        for k=j+1:numFailedChecks
            if ~strcmp(res.failedChecks(j).parameter,fieldname)
                j = k-1; break
            end
        end
    
        % list all infractions
        if j == i
            % single infraction
            fprintf("  - failed check for field '" + fieldname + "':\n");
            fprintf("    " + res.failedChecks(i).message + "\n");
        elseif j > i
            % multiple infractions
            fprintf("  - failed checks for field '" + fieldname + "':\n");
            for k=i:j
                fprintf("    " + res.failedChecks(k).message + "\n");
            end
        end
    end
end


% 3. redundant fields in struct for params/options

% read out redundant indices
redIdx = res.redIdxInUserList;

% simplify access to names of provided params / options
if strcmp(listname,'params')
    allfields = aux_readAllFields(params);
elseif strcmp(listname,'options')
    allfields = aux_readAllFields(options);
end

% if cond-default given, redIdx shorter than list -> expand
if length(allfields) > length(redIdx)
    redIdx = [redIdx; false(length(allfields) - length(redIdx),1)];
end

% init output text
redTxt = '''';

if (isa(sys,'hybridAutomaton') || isa(sys,'parallelHybridAutomaton')) ...
        && ~strcmp(func,'simulate')
    % simulateRandom -> simulateStandard (corresponding func in contDynamics)
    if isa(sys,'hybridAutomaton')
        % check contDynamics class from all locations
        for i=1:length(sys.location)
            classname = class(sys.location(i).contDynamics);
            if contains(func,'simulate')
                % for simulateRandom|simulate just contDynamics
                classname = 'contDynamics';
            end
            configfile{i} = ['config_' classname '_' func];
        end
    elseif isa(sys,'parallelHybridAutomaton')
        configfile = {};
        comps = zeros(length(sys.components),1);
        for i=1:length(sys.components)
            % check contDynamics class from all components/locations
            for j=1:length(sys.components(i).location)
                classname = class(sys.components(i).location(j).contDynamics);
                if contains(func,'simulate')
                    % for simulateNormal|simulate just contDynamics
                    classname = 'contDynamics';
                end
                configfile = [configfile; 'config_' classname '_' func];
            end
            comps(i) = j;
        end
        comps = cumsum(comps);
    end
    [configfile, cfgIdx] = unique(configfile);
    % read all configfile and add their list to hybrid list
    for i=1:length(configfile)
        if isa(sys,'hybridAutomaton')
            if strcmp(listname,'params')
                [contDynlist,~] = eval([configfile{i} '(sys.location(cfgIdx(i)).contDynamics,params,options)']);
            elseif strcmp(listname,'options')
                [~,contDynlist] = eval([configfile{i} '(sys.location(cfgIdx(i)).contDynamics,params,options)']);
            end
        elseif isa(sys,'parallelHybridAutomaton')
            cfgpHAcomp = find(cfgIdx(i)<=comps,1,'first');
            cfgpHAloc = cfgIdx(i) - comps(find(cfgIdx(i)>=comps,1,'last'));
            if isempty(cfgpHAloc)
                cfgpHAloc = cfgIdx(i);
            end
            if strcmp(listname,'params')
                [contDynlist,~] = eval([configfile{i} ...
                    '(sys.components(cfgpHAcomp).location(cfgpHAloc).contDynamics,params,options)']);
            elseif strcmp(listname,'options')
                [~,contDynlist] = eval([configfile{i} ...
                    '(sys.components(cfgpHAcomp).location(cfgpHAloc).contDynamics,params,options)']);
            end
        end
        % append new entries from contDynamics list
        idxNew = ~ismember({contDynlist.name}',{blueprint.name}');
        if any(idxNew)
            blueprint = [blueprint; contDynlist(idxNew)];
        end
    end
    % check if redIdx present in any of the added lists
    for i=1:length(redIdx)
        if redIdx(i)
            if any(ismember({blueprint.name}',allfields{i}))
                redIdx(i) = false;
            end
        end
    end
end

% first part of redundancy check: params / options which are not used
redTxt = [redTxt, strjoin(allfields(redIdx),''', ''')];

% second part of redundancy test: redundant conditional params / options
redIdx = ismember({blueprint.status}','redundant');
% check if redundant conditional params / options provided
for i=1:length(redIdx)
    if redIdx(i) && ~any(ismember(allfields,blueprint(i).name))
        redIdx(i) = false;
    end
end
% add comma for correct text output
if ~strcmp(redTxt,"'") && any(redIdx)
    redTxt = [redTxt, ', '];
end
% concatentate to first group
redTxt = [redTxt, strjoin({blueprint(redIdx).name},''', '''), ''''];

% print redundancies
if ~strcmp(redTxt,"''")
    if VALIDATEOPTIONS_ERRORS
        warning(['Redundant ' listname ': ' redTxt]);
    else
        headerPrinted = aux_printHeader(listname,headerPrinted);
        disp("  - redundant fields: " + redTxt);
    end
end

% additional message if mandatory parameters are missing or checks failed
if ~VALIDATEOPTIONS_ERRORS
    if any(res.mandMissing) || ~isempty(res.failedChecks)
        headerPrinted = aux_printHeader(listname,headerPrinted);
        fprintf("\n");
        disp("ALERT! Missing mandatory fields and/or failed " + ...
            "validation checks will likely cause run-time errors.");
        disp("Correct the wrong settings in order of appearance to resolve potential issues.");
    end

    % print footer (only if header printed, too)
    if headerPrinted
        disp("-*---------------------------------*-");
    end
end

end

function headerPrinted = aux_printHeader(listname,headerPrinted)
% helper function to print header at first instance of violations

if ~headerPrinted
    disp("-*---------------------------------*-");
    fprintf("Input argument validation (" + listname + "):\n");
    headerPrinted = true;
end

end

function allfields = aux_readAllFields(list)
% list      ... params struct / options struct

% auxiliary function to read out all fields from the params/options struct;
% we cannot only call 'fields(list)' since some parameters can also be
% struct, for which we therefore need extra manipulation. Finally, the
% output is a list with the format
%    'param1',
%    'param2.sub1',
%    'param2.sub2',
%    'param3', ...
% so that each field in a parameter-struct has its own (full) name

% read out all fields
allfields = fields(list);
% ... for each parameter that is itself a struct, 'allfields' now has only
% the name of that struct -> we want the full name though

% init index
idx = 0;

% increment index until all fields in 'allfields' are traversed
while idx < length(allfields)

    % increment index
    idx = idx + 1;

    % name of the parameter-struct
    structname = allfields{idx};

    % check if the entry in 'list' with the idx-th name in 'allfields' is
    % orginally a struct
    if isstruct(list.(structname))

        % list all fields of the parameter-struct
        structfields = fields(list.(structname));

        % shift that is required for all other parameters further down the
        % list to fit in the fields of the parameter-struct
        shift = length(structfields);

        % read out these other parameters
        temp = allfields(idx+1:end);

        % append names of fields of parameter-struct to 'allfields', using
        % the format 'param_i.sub_j'
        for j=1:shift
            allfields{idx-1+j} = [structname '.' structfields{j}];
        end
        
        % concatenate other parameter back onto the list of all parameters
        allfields = [allfields(1:idx-1+shift); temp];

        % correct index
        idx = idx + shift - 1;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
