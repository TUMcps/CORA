function [paramsList,optionsList] = config_contDynamics_simulateRandom(sys,params,options)
% config_contDynamics_simulateRandom - configuration file for validation
%    of model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_contDynamics_simulateRandom(sys,params,options)
%
% Inputs:
%    sys - contDynamics object
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
% Written:      03-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. init lists
initParamsOptionsLists();

% append entries to list of model parameters
add2params('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
add2params('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});
add2params('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
add2params('U','default',{@(val)any(ismember(getMembers('Usim'),class(val)))},{'memberUsim'});
add2params('u','default',{@isnumeric},{'isnumeric'});
add2params('tu','default',{@isvector,@isnumeric,@(val)withinTol(val(1),params.tStart),...
    @(val)length(val)==size(params.u,2)},...
    {'isvector','isnumeric','idx1eqtStart','equ'});
add2params('paramInt','mandatory',{@(val)length(val)==sys.nrOfParam,...
    @(val)isa(val,'interval') || (isvector(val) && isnumeric(val))},...
    {'eqparam','vectororinterval'},{@()isa(sys,'nonlinParamSys')});
add2params('y0guess','mandatory',{@(val)length(val)==sys.nrOfConstraints},...
    {'eqconstr'},{@()isa(sys,'nonlinDASys')});
add2params('W','default',{@(val)any(ismember(getMembers('W'),class(val))),...
    @(val)eq(dim(val),sys.dim)},{'memberW','eqsysdim'});
add2params('V','default',{@(val)any(ismember(getMembers('V'),class(val))),...
    @(val)eq(dim(val),sys.nrOfOutputs)},{'memberV','eqoutput'});

% append entries to list of algorithm parameters
add2options('type','default',{@ischar,@(val)any(ismember(getMembers('type'),val))},...
    {'ischar','membertype'});

% for all types = 'standard', 'gaussian', 'rrt'
add2options('points','default',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});
% only for types = 'standard' and 'gaussian'
add2options('nrConstInp','default',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,0),...
    @(val)c_nrConstInp(val,sys,params,options)},...
    {'isscalar','isnumeric','integer','gezero',''},...
    {@()any(strcmp(options.type,{'standard','gaussian'}))});
% only for type = 'standard'
add2options('fracVert','default',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'},{@()strcmp(options.type,'standard')});
add2options('fracInpVert','default',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'},{@()strcmp(options.type,'standard')});
% only for type = 'gaussian'
add2options('p_conf','default',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'},{@()strcmp(options.type,'gaussian')});
% only for type = 'rrt'
add2options('vertSamp','mandatory',{@isscalar,@islogical},{'isscalar','islogical'},...
    {@()strcmp(options.type,'rrt')});
add2options('stretchFac','mandatory',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},{@()strcmp(options.type,'rrt')});
add2options('R','mandatory',{@(val)isa(val,'reachSet')},{'isareachSet'},...
    {@()strcmp(options.type,'rrt')});


% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------
