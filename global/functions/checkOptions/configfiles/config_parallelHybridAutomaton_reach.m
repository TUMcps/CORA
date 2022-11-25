function [paramsList,optionsList] = config_parallelHybridAutomaton_reach(sys,params,options)
% config_parallelhybridAutomaton_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_parallelhybridAutomaton_reach
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
% Last revision:---

%------------- BEGIN CODE --------------

% 1. init lists
initParamsOptionsLists();

% append entries to list of model parameters
add2params('tStart','default',{@isscalar,@isnumeric,@(val)ge(val,0)},{'isscalar','isnumeric','gezero'});
add2params('tFinal','mandatory',{@isscalar,@isnumeric,@(val)val>0},{'isscalar','isnumeric','gezero'});
add2params('startLoc','mandatory',{@(val)c_pHA_startLoc(val,sys,params)},{''});
add2params('finalLoc','default',{@(val)c_pHA_finalLoc(val,sys,params)},{''});
add2params('R0','mandatory',{@(val)c_pHA_R0(val,sys,params)},{''});
add2params('inputCompMap','default',...
    {@isnumeric,@(val)max(val)<=length(sys.components),@(val)min(val)>=1,...
    @(val) all(size(val)==[sys.numInputs,1]) || all(size(val)==[1,sys.numInputs])},...
    {'isnumeric','lecomp','vectorgeone','???'});
add2params('U','default',{@(val)c_pHA_U(val,sys,params)},{''});


% append entries to list of algorithm parameters
add2options('timeStep','mandatory',{@(val)c_HA_timeStep(val,sys,options)},{''});
add2options('guardIntersect','mandatory',{@ischar,@(val)any(ismember(getMembers('guardIntersect'),val))},...
    {'ischar','memberguardIntersect'});
add2options('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
add2options('enclose','mandatory',{@iscell,@(val)any(ismember(getMembers('enclose'),val))},...
    {'iscell','memberenclose'},{@()any(ismember(getMembers('guardIntersect4enclose'),options.guardIntersect))});
add2options('guardOrder','mandatory',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},{@()any(ismember(getMembers('guardIntersect4guardOrder'),options.guardIntersect))});
add2options('intersectInvariant','optional',{@isscalar,@islogical},{'isscalar','islogical'});

% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------

