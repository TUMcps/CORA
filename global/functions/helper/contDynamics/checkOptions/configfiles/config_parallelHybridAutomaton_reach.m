function [paramsList,optionsList] = config_parallelHybridAutomaton_reach(sys,params,options)
% config_parallelHybridAutomaton_reach - configuration file for validation
%    of model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = ...
%       config_parallelHybridAutomaton_reach(sys,params,options)
%
% Inputs:
%    sys - parallelHybridAutomaton object
%    params - user-defined model parameters
%    options - user-defined algorithm parameters
%
% Outputs:
%    paramsList - list of model parameters
%    optionsList - list of algorithm parameters
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      04-February-2021
% Last update:  ---
% Last revision:19-June-2023 (MW, structs, remove global variables)

%------------- BEGIN CODE --------------

% list of model parameters
paramsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of model parameters
paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@isnumeric,@(val)ge(val,0)},{'isscalar','isnumeric','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@isnumeric,@(val)val>0},{'isscalar','isnumeric','gezero'});
paramsList(end+1,1) = add2list('startLoc','mandatory',{@(val)c_pHA_startLoc(val,sys,params)},{''});
paramsList(end+1,1) = add2list('finalLoc','default',{@(val)c_pHA_finalLoc(val,sys,params)},{''});
paramsList(end+1,1) = add2list('R0','mandatory',{@(val)c_pHA_R0(val,sys,params)},{''});
paramsList(end+1,1) = add2list('inputCompMap','default',...
    {@isnumeric,@(val)max(val)<=length(sys.components),@(val)min(val)>=1,...
    @(val) all(size(val)==[sys.nrOfInputs,1]) || all(size(val)==[1,sys.nrOfInputs])},...
    {'isnumeric','lecomp','vectorgeone','???'});
paramsList(end+1,1) = add2list('U','default',{@(val)c_pHA_U(val,sys,params)},{''});
paramsList(end+1,1) = add2list('u','default',...
    {@iscell,@(val)all(size(val) == [length(sys.components),1])},{'iscell','???'});

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of algorithm parameters
optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('timeStep','mandatory',{@(val)c_HA_timeStep(val,sys,options)},{''});
optionsList(end+1,1) = add2list('guardIntersect','mandatory',{@ischar,@(val)any(ismember(getMembers('guardIntersect'),val))},...
    {'ischar','memberguardIntersect'});
optionsList(end+1,1) = add2list('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
optionsList(end+1,1) = add2list('enclose','mandatory',{@iscell,@(val)any(ismember(getMembers('enclose'),val))},...
    {'iscell','memberenclose'},{@()any(ismember(getMembers('guardIntersect4enclose'),options.guardIntersect))});
optionsList(end+1,1) = add2list('guardOrder','mandatory',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},{@()any(ismember(getMembers('guardIntersect4guardOrder'),options.guardIntersect))});
optionsList(end+1,1) = add2list('intersectInvariant','default',{@isscalar,@islogical},{'isscalar','islogical'});

end

%------------- END OF CODE --------------
