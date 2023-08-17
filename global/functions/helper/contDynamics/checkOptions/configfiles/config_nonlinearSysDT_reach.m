function [paramsList,optionsList] = config_nonlinearSysDT_reach(sys,params,options)
% config_nonlinearSysDT_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_nonlinearSysDT_reach(sys,params,options)
%
% Inputs:
%    sys - nonlinearSysDT object
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
paramsList(end+1,1) = add2list('U','default',{@(val)any(ismember(getMembers('U'),class(val)))},{'memberU'});
paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});
paramsList(end+1,1) = add2list('u','default',{@isnumeric},{'isnumeric',''});

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of algorithm parameters
optionsList(end+1,1) = add2list('alg','default',{@ischar,@(val)any(ismember(getMembers('alg4DT'),val))},{'ischar','memberalg4DT'});
optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('zonotopeOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'},...
    {@()~contains(options.alg,'adaptive')});
optionsList(end+1,1) = add2list('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
optionsList(end+1,1) = add2list('tensorOrder','mandatory',...
    {@isscalar,@(val)mod(val,1)==0,@(val)any(val==[2,3])},{'isscalar','integer','2or3'},...
    {@()~contains(options.alg,'adaptive')});
optionsList(end+1,1) = add2list('errorOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'},...
    {@()~contains(options.alg,'adaptive')});
optionsList(end+1,1) = add2list('tensorOrderOutput','default',...
    {@isscalar,@(val)mod(val,1)==0,@(val)any(val==[2,3])},...
    {'isscalar','integer','2or3'});

optionsList(end+1,1) = add2list('compOutputSet','default',{@isscalar,@islogical},...
    {'isscalar','islogical'});

% polyZono
optionsList(end+1,1) = add2list('polyZono.maxDepGenOrder','default',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},...
    {@()isa(params.R0,'polyZonotope')});
optionsList(end+1,1) = add2list('polyZono.maxPolyZonoRatio','default',{@isscalar,@isnumeric,@(val)val>0},...
    {'isscalar','isnumeric','gezero'},...
    {@()isa(params.R0,'polyZonotope')});
optionsList(end+1,1) = add2list('polyZono.restructureTechnique','default',{@ischar,...
    @(val)any(ismember(getMembers('restructureTechnique'),val))},...
    {'ischar','memberrestructureTechnique'},...
    {@()isa(params.R0,'polyZonotope')});

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