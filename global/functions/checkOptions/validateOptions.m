function options = validateOptions(sys,func,params,options)
% validateOptions - validates model parameters and algorithm parameters
%    for given function
%
% Syntax:
%    options = validateOptions(sys,func,params,options)
%
% Inputs:
%    sys - object of contDynamics class
%    func - function for which parameters should be validated
%    params - model parameters (e.g., initial set)
%    options - algorithm parameters (e.g., time step size)
%
% Outputs:
%    options - unified parameters for further computations
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      25-January-2021
% Last update:  26-January-2021
% Last revision:---

%------------- BEGIN CODE --------------

% see if configuration file for given class + func exists
currfolder = mfilename('fullpath');
currfolder = currfolder(1:end-length(mfilename));
classname = class(sys);
% rename classname to contDynamics for any simulate or observe function
if isa(sys,'contDynamics') && any(contains(func,{'simulate','observe'}))
    classname = 'contDynamics';
end
configfile = ['config_' classname '_' func];
if ~isfile([currfolder filesep 'configfiles' filesep configfile '.m'])
    error("No configuration file for desired validation file");
end

% initialize codex of error messages (global variable)
initErrorCodex; % access with: 'global codex;'

% check params
params = checkFullList(sys,func,params,options,configfile,'params');

% check options
options = checkFullList(sys,func,params,options,configfile,'options');

% post-processing: set params / options
[params, options] = postProcessing(sys,func,params,options);

% unify params and options
options = params2options(params,options);

% clear codex of error messages (global variable)
clear codex;


end


% Auxiliary Functions -----------------------------------------------------

function list = checkFullList(sys,func,params,options,configfile,listname)
% list          ... params struct / options struct
% blueprint     ... paramsList / optionsList (from config file)
% sys           ... object of some contDynamics class
% func          ... name of checked function, i.e. 'reach'
% configfile    ... name of configuration file
% listname      ... 'params', 'options', for error messages

% read configuration file
if strcmp(listname,'params')
    [blueprint,~] = eval([configfile '(sys,params,options)']);
    list = params;
elseif strcmp(listname,'options')
    [~,blueprint] = eval([configfile '(sys,params,options)']);
    list = options;
end

% set non-cond defaults
list = setMissingDefaults(list,blueprint,sys,listname);

% update configuration file
if strcmp(listname,'params')
    [blueprint,~] = eval([configfile '(sys,list,options)']);
elseif strcmp(listname,'options')
    [~,blueprint] = eval([configfile '(sys,params,list)']);
end

% read all fields (special read out when field is a struct)
allfields = readAllFields(list);

% redundant indices
redIdx = false(length(allfields),1);

% conditional indices
condIdx = false(length(allfields),1);

% go through all fields, see which are conditional or redundant
for i=1:length(allfields)

    % find idx of current field in blueprint
    idx = find(ismember(blueprint.name,allfields{i}),1,'first');

    if isempty(idx)
        % field not in blueprint -> param / option is redundant
        redIdx(i) = true;
    elseif ~isempty(blueprint.condfunc{idx})
        % field not in blueprint -> param / option is redundant
        condIdx(i) = true;
    end
end

% skip conditional and redundant indices
skipIdx = condIdx | redIdx;

% check non-conditional and non-redundant params / options
validateList(list,blueprint,sys,listname,skipIdx);

% check conditions for conditional params / options
% ... remove 'cond-' from status
[blueprint, condBPIdx, condIdx] = checkCond(list,blueprint);

if any(condBPIdx)
    % set conditional defaults (now also listed in blueprint as 'default')
    list = setMissingDefaults(list,blueprint,sys,listname);
    
    % another update of blueprint would be required if any of the
    % conditional defaults were set automatically to their default value
    % and they would be in a check function of another conditional param / option  
    
    % unify condIdx
    condIdx = [condIdx; true(length(readAllFields(list))-length(condIdx),1)];
    
    % remaining: check conditional params / options (also cond-default...)
    validateList(list,blueprint,sys,listname,~condIdx);
end

% print redundancies
if strcmp(listname,'params')
    redundantFields(list,options,blueprint,sys,func,listname,redIdx);
elseif strcmp(listname,'options')
    redundantFields(list,params,blueprint,sys,func,listname,redIdx);
end


end

function list = setMissingDefaults(list,blueprint,sys,listname)
% list      ... params struct / options struct
% blueprint ... paramsList / optionsList (from config file)
% sys       ... object of some hybrid or contDynamics class
% listname  ... 'params' or 'options'

% only checked whether default params / options are there or not
% if not -> add that param / option with default value

% indices with default values (skip where condfunc given)
idxDefault = ismember(blueprint.status,'default') & ...
    cellfun(@(x)isempty(x),blueprint.condfunc);
% default fields
defFields = blueprint.name(idxDefault);
% names of user-defined fields
fieldnames = readAllFields(list);
% loop over default params / options
for i=1:length(defFields)
    % check if default param / option given
    if ~any(ismember(fieldnames,defFields{i}))
        % assign default value by reading from list of default values
        if ~contains(defFields{i},'.')
            % standard param / option
            list.(defFields{i}) = getDefaultValue(defFields{i},sys,listname);
        else
            % param / option is a struct
            dotIdx = strfind(defFields{i},'.');
            firstname = defFields{i}(1:dotIdx-1);
            secondname = defFields{i}(dotIdx+1:end);
            list.(firstname).(secondname) = getDefaultValue(defFields{i},sys,listname);
        end
    end
end

end

function validateList(list,blueprint,sys,listname,skipIdx)
% list      ... params struct / options struct
% blueprint ... paramsList / optionsList (from config file)
% sys       ... object of some hybrid or contDynamics class
% listname  ... 'params', 'options', for error messages
% skipIdx   ... indices in list which are to be skipped
%               (redundant and conditional params / options)

% read fields correctly from list
allfields = readAllFields(list);

% containers
checkSuccessful = false(length(allfields),1);
checkSuccessful(skipIdx) = true;
sanitycounter = 0;

% check mandatory status of blueprint
for i=1:length(blueprint.name)
    if strcmp(blueprint.status{i},'mandatory') && isempty(blueprint.condfunc{i})
        if ~isfield(list,blueprint.name{i})
            % mandatory option missing -> print error
            error(sprintf(...
                'Error in %s check for %s object:\n  %s.%s is missing.',...
                listname,class(sys),listname,blueprint.name{i}));
        end
    end
end

while ~all(checkSuccessful)
    
    % go through all parameters (also defaults) in correct order
    for i=1:length(blueprint.name)
        
        % find idx of current field in blueprint
        idx = find(ismember(allfields,blueprint.name{i}),1,'first');
        
        % skip if already successful
        if ~isempty(idx) && ~checkSuccessful(idx)
            
            % special handling for struct
            if contains(blueprint.name{i},'.')
                dotIdx = strfind(blueprint.name{i},'.');
                firstname = blueprint.name{i}(1:dotIdx-1);
                secondname = blueprint.name{i}(dotIdx+1:end);
            end
            
            % check validity of assigned value
            for j=1:length(blueprint.checkfuncs{i})
                % current method
%                 if ~contains(blueprint.name{i},'.') % standard
%                     if ~blueprint.checkfuncs{i}{j}(list.(blueprint.name{i}))
%                         error(printOptionOutOfRange(sys,blueprint.name{i},listname));
%                     end
%                 else % special handling for struct
%                     if ~blueprint.checkfuncs{i}{j}(list.(firstname).(secondname))
%                         error(printOptionOutOfRange(sys,blueprint.name{i},listname));
%                     end
%                 end
                
                % future handling of error messages
                errmsgid = blueprint.errmsgs{i}{j};
                if ~contains(blueprint.name{i},'.') % standard
                    if isempty(errmsgid) % c_* validation functions return also msg
                        [rescheck,msg] = blueprint.checkfuncs{i}{j}(list.(blueprint.name{i}));
                    else 
                        rescheck = blueprint.checkfuncs{i}{j}(list.(blueprint.name{i}));
                        msg = getErrorMessage(errmsgid);
                    end
                else % special handling for struct
                    if isempty(errmsgid) % c_* validation functions return also msg
                        [rescheck,msg] = blueprint.checkfuncs{i}{j}(list.(firstname).(secondname));
                    else 
                        rescheck = blueprint.checkfuncs{i}{j}(list.(firstname).(secondname));
                        msg = getErrorMessage(errmsgid);
                    end
                end
                % error message if validation failed
                if ~rescheck
                    error(sprintf([listname '.' blueprint.name{i} ' ' msg '.']));
                end
            end
            % all checks for parameter have been successful
            checkSuccessful(idx) = true;
            
        end
    end
    
    sanitycounter = sanitycounter + 1;
    if sanitycounter == 5 % prevent infinite loop
        error("Bug in definition of validation function. Report to devs.");
    end

end

end

function [blueprint,condBPIdx,condIdx] = checkCond(list,blueprint)
% list      ... params struct / options struct
% blueprint ... paramsList / optionsList (from config file)

condBPIdx = false(length(blueprint.name),1);
for i=1:length(blueprint.name)
    if ~isempty(blueprint.condfunc{i})
        condBPIdx(i) = true;
        % check condition
        if ~blueprint.condfunc{i}()
            % set field 'status' to 'redundant'
            blueprint.status{i} = 'redundant';
            condBPIdx(i) = false;
        end
        % remove condition from validation functions
        blueprint.condfunc{i} = {};
    end
end

% find condBPIdx entries in list
allfields = readAllFields(list);
condIdx = false(length(allfields),1);
for i=1:length(condBPIdx)
    if condBPIdx(i)
        idx = find(ismember(allfields,blueprint.name{i}),1,'first');
        condIdx(idx) = true;
    end
end

end

function redundantFields(list,otherlist,blueprint,sys,func,listname,redIdx)
% list      ... params struct / options struct
% otherlist ... complementary: params struct / options struct
% blueprint ... paramsList / optionsList (from config file)
% sys       ... object of some hybrid or contDynamics class
% func      ... name of checked function, i.e. 'reach'
% listname  ... 'params' or 'options'
% redIdx    ... redundant indices from original list

% simplify access to names of provided params / options
allfields = readAllFields(list);

% if cond-default given, redIdx shorter than list -> expand
if length(allfields) > length(redIdx)
    redIdx = [redIdx; false(length(allfields) - length(redIdx),1)];
end

% init output text
redTxt = '';

if (isa(sys,'hybridAutomaton') || isa(sys,'parallelHybridAutomaton')) ...
        && ~strcmp(func,'simulate')
    % simulateRandom -> simulateStandard (corresponding func in contDynamics)
    if strcmp(func,'simulateRandom')
        func = 'simulateStandard';
    end
    if isa(sys,'hybridAutomaton')
        % check contDynamics class from all locations
        for i=1:length(sys.location)
            classname = class(sys.location{i}.contDynamics);
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
            for j=1:length(sys.components{i}.location)
                classname = class(sys.components{i}.location{j}.contDynamics);
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
                [contDynlist,~] = eval([configfile{i} '(sys.location{cfgIdx(i)}.contDynamics,list,otherlist)']);
            elseif strcmp(listname,'options')
                [~,contDynlist] = eval([configfile{i} '(sys.location{cfgIdx(i)}.contDynamics,otherlist,list)']);
            end
        elseif isa(sys,'parallelHybridAutomaton')
            cfgpHAcomp = find(cfgIdx(i)<=comps,1,'first');
            cfgpHAloc = cfgIdx(i) - comps(find(cfgIdx(i)>=comps,1,'last'));
            if isempty(cfgpHAloc); cfgpHAloc = cfgIdx(i); end
            if strcmp(listname,'params')
                [contDynlist,~] = eval([configfile{i} ...
                    '(sys.components{cfgpHAcomp}.location{cfgpHAloc}.contDynamics,list,otherlist)']);
            elseif strcmp(listname,'options')
                [~,contDynlist] = eval([configfile{i} ...
                    '(sys.components{cfgpHAcomp}.location{cfgpHAloc}.contDynamics,otherlist,list)']);
            end
        end
        blueprint.name = [blueprint.name; contDynlist.name];
        blueprint.status = [blueprint.status; contDynlist.status];
        blueprint.checkfuncs = [blueprint.checkfuncs; contDynlist.checkfuncs];
    end
    % check if redIdx present in any of the added lists
    for i=1:length(redIdx)
        if redIdx(i)
            if any(ismember(blueprint.name,allfields{i}))
                redIdx(i) = false;
            end
        end
    end
end

% first part of redundancy check: params / options which are not used
redTxt = [redTxt, strjoin(allfields(redIdx),', ')];

% second part of redundancy test: redundant conditional params / options
redIdx = ismember(blueprint.status,'redundant');
% check if redundant conditional params / options provided
for i=1:length(redIdx)
    if redIdx(i) && ~any(ismember(allfields,blueprint.name(i)))
        redIdx(i) = false;
    end
end
% add comma for correct text output
if ~isempty(redTxt) && any(redIdx)
    redTxt = [redTxt, ', '];
end
% concatentate to first group
redTxt = [redTxt, strjoin(blueprint.name(redIdx),', ')];

% print redundancies
if ~isempty(redTxt)
    warning(['Redundant ' listname ': ' redTxt]);
end

end

function allfields = readAllFields(list)
% list      ... params struct / options struct

allfields = fields(list);
idx = 0;
while idx < length(allfields)
    idx = idx + 1;
    if isstruct(list.(allfields{idx}))
        structfields = fields(list.(allfields{idx}));
        shift = length(structfields);
        temp = allfields(idx+1:end);
        structname = allfields{idx};
        for j=1:shift
            allfields{idx-1+j} = [structname '.' structfields{j}];
        end
        allfields = [allfields(1:idx-1+shift); temp];
        idx = idx + shift - 1;
    end
    if idx == 100
        error("Bug!"); % remove...
    end
end

end

%------------- END OF CODE --------------

