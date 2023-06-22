function [paramsList,optionsList] = config_nonlinearSys_reachInnerProjection(sys,params,options)
% config_nonlinearSys_reachInnerProjection - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_nonlinearSys_reachInnerProjection(sys,params,options)
%
% Inputs:
%    sys - nonlinearSys object
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
paramsList(end+1,1) = add2list('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of algorithm parameters
optionsList(end+1,1) = add2list('algInner','mandatory',{@ischar,@(val)any(ismember(getMembers('algInner'),val))},...
    {'ischar','memberalgInner'});

optionsList(end+1,1) = add2list('timeStep','mandatory',{@isscalar,@(val)val>0,...
    @(val)abs(params.tFinal/val - round(params.tFinal/val))<1e-9,...
	@(val)c_inputTraj(val,sys,params,options)},{'isscalar','gezero','intsteps',''});
optionsList(end+1,1) = add2list('taylorOrder','mandatory',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});
optionsList(end+1,1) = add2list('taylmOrder','mandatory',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});

optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});

% polyZono
optionsList(end+1,1) = add2list('polyZono.maxDepGenOrder','default',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'});
optionsList(end+1,1) = add2list('polyZono.maxPolyZonoRatio','default',{@isscalar,@isnumeric,@(val)val>0},...
    {'isscalar','isnumeric','gezero'});
optionsList(end+1,1) = add2list('polyZono.restructureTechnique','default',{@ischar,...
    @(val)any(ismember(getMembers('restructureTechnique'),val))},...
    {'ischar','memberrestructureTechnique'});

% lagrangeRem
optionsList(end+1,1) = add2list('lagrangeRem.simplify','default',{@ischar,...
    @(val)any(ismember(getMembers('lagrangeRem.simplify'),val))},...
    {'ischar','memberlagrangeRem.simplify'});
optionsList(end+1,1) = add2list('lagrangeRem.method','default',{@ischar,...
    @(val)any(ismember(getMembers('lagrangeRem.method'),val))},...
    {'ischar','memberlagrangeRem.method'});
optionsList(end+1,1) = add2list('lagrangeRem.tensorParallel','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('lagrangeRem.replacements','optional',{@(val)isa(val,'function_handle')},{'isafunction_handle'});
optionsList(end+1,1) = add2list('lagrangeRem.zooMethods','mandatory',...
    {@ischar,@(val)any(ismember(getMembers('lagrangeRem.zooMethods'),val))},...
    {'ischar','memberlagrangeRem.zooMethods'},...
    {@()strcmp(options.lagrangeRem.method,'zoo')});
optionsList(end+1,1) = add2list('lagrangeRem.optMethod','default',...
    {@ischar,@(val)any(ismember(getMembers('lagrangeRem.optMethod'),val))},...
    {'ischar','memberlagrangeRem.optMethod'},...
    {@()strcmp(options.lagrangeRem.method,'taylorModel')});
optionsList(end+1,1) = add2list('lagrangeRem.maxOrder','optional',...
    {@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'},...
    {@()any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}))});
optionsList(end+1,1) = add2list('lagrangeRem.tolerance','optional',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},...
    {@()any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}))});
optionsList(end+1,1) = add2list('lagrangeRem.eps','optional',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},...
    {@()any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}))});

end

%------------- END OF CODE --------------

