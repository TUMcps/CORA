function [paramsList,optionsList] = config_hybridAutomaton_reach(sys,params,options)
% config_hybridAutomaton_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_hybridAutomaton_reach(sys,params,options)
%
% Inputs:
%    sys - linParamSys object
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
% Written:      27-January-2021
% Last update:  04-February-2021
% Last revision:19-June-2023 (MW, structs, remove global variables)

%------------- BEGIN CODE --------------

% list of model parameters
paramsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of model parameters
paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@isnumeric,@(val)ge(val,0)},{'isscalar','isnumeric','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@isnumeric,@(val)val>0},{'isscalar','isnumeric','gezero'});
paramsList(end+1,1) = add2list('startLoc','mandatory',{@isscalar,@isnumeric,...
    @(val)mod(val,1)==0,@(val)ge(val,1),@(val)le(val,length(sys.location))},...
    {'isscalar','isnumeric','integer','geone','leloc'});
paramsList(end+1,1) = add2list('finalLoc','default',{@isscalar,@isnumeric,...
    @(val)mod(val,1)==0,@(val)ge(val,0),@(val)le(val,length(sys.location)+1)},...
    {'isscalar','isnumeric','integer','geone','lelocplus1'});
paramsList(end+1,1) = add2list('R0','mandatory',{@(val)isa(val,'contSet'),...
    @(val)dim(val)==sys.location(params.startLoc).contDynamics.dim},{'isacontSet','eqsysdim'});
paramsList(end+1,1) = add2list('U','default',{@(val) (~iscell(val) && isa(val,'contSet') ) || ...
    ( iscell(val) && (all(size(val) == [length(sys.location),1]) || ...
        all(size(val) == [1,length(sys.location)])) )},{'isacontSet','???'});

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of algorithm parameters
optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('linAlg','default',{@ischar,@(val)any(ismember(getMembers('linAlg4HA'),val))},...
    {'ischar','memberlinAlg4HA'});
optionsList(end+1,1) = add2list('timeStep','mandatory',{@(val)c_HA_timeStep(val,sys,options)},...
    {''},{@()~strcmp(options.linAlg,'adaptive')});
optionsList(end+1,1) = add2list('error','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},{@()strcmp(options.linAlg,'adaptive')});
optionsList(end+1,1) = add2list('guardIntersect','mandatory',{@ischar,@(val)any(ismember(getMembers('guardIntersect'),val))},...
    {'ischar','memberguardIntersect'});
optionsList(end+1,1) = add2list('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
optionsList(end+1,1) = add2list('enclose','mandatory',{@iscell,@(val)any(ismember(getMembers('enclose'),val))},...
    {'iscell','memberenclose'},{@()any(ismember(getMembers('guardIntersect4enclose'),options.guardIntersect))});
optionsList(end+1,1) = add2list('guardOrder','mandatory',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},...
    {@()any(ismember(getMembers('guardIntersect4guardOrder'),options.guardIntersect))});
optionsList(end+1,1) = add2list('intersectInvariant','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('compTimePoint','default',{@isscalar,@islogical},{'isscalar','islogical'});

end

%------------- END OF CODE --------------

