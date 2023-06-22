function [paramsList,optionsList] = config_parallelHybridAutomaton_simulate(sys,params,options)
% config_parallelHybridAutomaton_simulate - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_parallelHybridAutomaton_simulate(sys,params,options)
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
paramsList(end+1,1) = add2list('x0','mandatory',{@(val)all(size(val)==[sys.dim,1])},{'eqsysdim'});
paramsList(end+1,1) = add2list('inputCompMap','default',...
    {@isnumeric,@(val)max(val)<=length(sys.components),@(val)min(val)>=1,...
    @(val) all(size(val)==[sys.nrOfInputs,1]) || all(size(val)==[1,sys.nrOfInputs])},...
    {'isnumeric','lecomp','vectorgeone','???'});
paramsList(end+1,1) = add2list('u','default',{@(val)c_pHA_sim_u(val,sys,params)},{''});

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

end

%------------- END OF CODE --------------


