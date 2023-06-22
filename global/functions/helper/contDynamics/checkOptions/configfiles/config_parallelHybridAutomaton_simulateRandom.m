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

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of algorithm parameters
optionsList(end+1,1) = add2list('points','default',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});
optionsList(end+1,1) = add2list('fracVert','default',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'});
optionsList(end+1,1) = add2list('fracInpVert','default',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'});
optionsList(end+1,1) = add2list('nrConstInp','default',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});

end

%------------- END OF CODE --------------


