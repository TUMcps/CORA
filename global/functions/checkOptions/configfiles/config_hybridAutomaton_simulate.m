function [paramsList,optionsList] = config_hybridAutomaton_simulate(sys,params,options)
% config_hybridAutomaton_simulate - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_hybridAutomaton_simulate(sys,params,options)
%
% Inputs:
%    sys - hybridAutomaton object
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
add2params('startLoc','mandatory',{@isscalar,@isnumeric,...
    @(val)mod(val,1)==0,@(val)ge(val,1),@(val)le(val,length(sys.location))},...
    {'isscalar','isnumeric','integer','geone','leloc'});
add2params('finalLoc','default',{@isscalar,@isnumeric,...
    @(val)mod(val,1)==0,@(val)ge(val,0),@(val)le(val,length(sys.location))},...
    {'isscalar','isnumeric','integer','geone','leloc'});
add2params('x0','mandatory',{@isvector,@(val)length(val)==sys.location{params.startLoc}.contDynamics.dim},...
    {'isvector','eqsysdim'});
add2params('u','default',{@(val)c_HA_sim_u(val,sys,params)},{''});

% append entries to list of algorithm parameters
% ---

% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------


