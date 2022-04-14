function [paramsList,optionsList] = config_parallelHybridAutomaton_simulateRandom(sys,params,options)
% config_parallelHybridAutomaton_simulateRandom - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_parallelHybridAutomaton_simulateRandom(sys,params,options)
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
add2options('points','mandatory',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});
add2options('fracVert','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'});
add2options('fracInpVert','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'});
add2options('inpChanges','mandatory',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,0)},...
    {'isscalar','isnumeric','integer','gezero'});

% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------


