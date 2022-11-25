function [paramsList,optionsList] = config_hybridAutomaton_reach(sys,params,options)
% config_hybridAutomaton_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_hybridAutomaton_reach
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
% Last revision:---

%------------- BEGIN CODE --------------

% 1. init lists
initParamsOptionsLists();

% append entries to list of model parameters
add2params('tStart','default',{@isscalar,@isnumeric,@(val)ge(val,0)},{'isscalar','isnumeric','gezero'});
add2params('tFinal','mandatory',{@isscalar,@isnumeric,@(val)val>0},{'isscalar','isnumeric','gezero'});
add2params('startLoc','mandatory',{@isscalar,@isnumeric,...
    @(val)mod(val,1)==0,@(val)ge(val,1),@(val)le(val,length(sys.location))},...
    {'isscalar','isnumeric','integer','geone','leloc'});
add2params('finalLoc','default',{@isscalar,@isnumeric,...
    @(val)mod(val,1)==0,@(val)ge(val,0),@(val)le(val,length(sys.location)+1)},...
    {'isscalar','isnumeric','integer','geone','lelocplus1'});
add2params('R0','mandatory',{@(val)isa(val,'contSet'),...
    @(val)dim(val)==sys.location{params.startLoc}.contDynamics.dim},{'isacontSet','eqsysdim'});
add2params('U','default',{@(val) (~iscell(val) && isa(val,'contSet') ) || ...
    ( iscell(val) && (all(size(val) == [length(sys.location),1]) || ...
        all(size(val) == [1,length(sys.location)])) )},{'isacontSet','???'});

% append entries to list of algorithm parameters
add2options('linAlg','default',{@ischar,@(val)any(ismember(getMembers('linAlg4HA'),val))},...
    {'ischar','memberlinAlg4HA'});
add2options('timeStep','mandatory',{@(val)c_HA_timeStep(val,sys,options)},...
    {''},{@()~strcmp(options.linAlg,'adaptive')});
add2options('error','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},{@()strcmp(options.linAlg,'adaptive')});
add2options('guardIntersect','mandatory',{@ischar,@(val)any(ismember(getMembers('guardIntersect'),val))},...
    {'ischar','memberguardIntersect'});
add2options('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
add2options('enclose','mandatory',{@iscell,@(val)any(ismember(getMembers('enclose'),val))},...
    {'iscell','memberenclose'},{@()any(ismember(getMembers('guardIntersect4enclose'),options.guardIntersect))});
add2options('guardOrder','mandatory',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},...
    {@()any(ismember(getMembers('guardIntersect4guardOrder'),options.guardIntersect))});
add2options('intersectInvariant','optional',{@isscalar,@islogical},{'isscalar','islogical'});
add2options('compTimePoint','default',{@isscalar,@islogical},{'isscalar','islogical'});

% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------

