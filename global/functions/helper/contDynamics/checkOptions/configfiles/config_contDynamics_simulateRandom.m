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
% Last revision:19-June-2023 (MW, structs, remove global variables)

%------------- BEGIN CODE --------------

% list of model parameters
paramsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of model parameters
paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});
paramsList(end+1,1) = add2list('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
paramsList(end+1,1) = add2list('U','default',{@(val)any(ismember(getMembers('Usim'),class(val)))},{'memberUsim'});
paramsList(end+1,1) = add2list('u','default',{@isnumeric},{'isnumeric'});
paramsList(end+1,1) = add2list('tu','default',{@isvector,@isnumeric,@(val)withinTol(val(1),params.tStart),...
    @(val)length(val)==size(params.u,2)},...
    {'isvector','isnumeric','idx1eqtStart','equ'});
paramsList(end+1,1) = add2list('paramInt','mandatory',{@(val)length(val)==sys.nrOfParam,...
    @(val)isa(val,'interval') || (isvector(val) && isnumeric(val))},...
    {'eqparam','vectororinterval'},{@()isa(sys,'nonlinParamSys')});
paramsList(end+1,1) = add2list('y0guess','mandatory',{@(val)length(val)==sys.nrOfConstraints},...
    {'eqconstr'},{@()isa(sys,'nonlinDASys')});
paramsList(end+1,1) = add2list('W','default',{@(val)any(ismember(getMembers('W'),class(val))),...
    @(val)eq(dim(val),sys.dim)},{'memberW','eqsysdim'});
paramsList(end+1,1) = add2list('V','default',{@(val)any(ismember(getMembers('V'),class(val))),...
    @(val)eq(dim(val),sys.nrOfOutputs)},{'memberV','eqoutput'});


% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

optionsList(end+1,1) = add2list('type','default',{@ischar,@(val)ismember(getMembers('type'),val)},...
    {'ischar','membertype'});
% for all types = 'standard', 'gaussian', 'rrt'
optionsList(end+1,1) = add2list('points','default',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});
% only for types = 'standard' and 'gaussian'
optionsList(end+1,1) = add2list('nrConstInp','default',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,0),...
    @(val)c_nrConstInp(val,sys,params,options)},...
    {'isscalar','isnumeric','integer','gezero',''},...
    {@()any(strcmp(options.type,{'standard','gaussian'}))});
% only for type = 'standard'
optionsList(end+1,1) = add2list('fracVert','default',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'},{@()strcmp(options.type,'standard')});
optionsList(end+1,1) = add2list('fracInpVert','default',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'},{@()strcmp(options.type,'standard')});
% only for type = 'gaussian'
optionsList(end+1,1) = add2list('p_conf','default',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'},{@()strcmp(options.type,'gaussian')});
% only for type = 'rrt'
optionsList(end+1,1) = add2list('vertSamp','mandatory',{@isscalar,@islogical},{'isscalar','islogical'},...
    {@()strcmp(options.type,'rrt')});
optionsList(end+1,1) = add2list('stretchFac','mandatory',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},{@()strcmp(options.type,'rrt')});
optionsList(end+1,1) = add2list('R','mandatory',{@(val)isa(val,'reachSet')},{'isareachSet'},...
    {@()strcmp(options.type,'rrt')});

end

%------------- END OF CODE --------------
