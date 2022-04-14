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
% Last revision:---

%------------- BEGIN CODE --------------

% 1. init lists
initParamsOptionsLists();

% append entries to list of model parameters
add2params('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
add2params('U','default',{@(val)any(ismember(getMembers('U'),class(val)))},{'memberU'});
add2params('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
add2params('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});
add2params('u','default',{@isnumeric},{'isnumeric',''});

% append entries to list of algorithm parameters
add2options('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
add2options('zonotopeOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});
add2options('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
add2options('tensorOrder','mandatory',...
    {@isscalar,@(val)mod(val,1)==0,@(val)any(val==[2,3])},{'isscalar','integer','2or3'});
add2options('errorOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});

% polyZono
add2options('polyZono.maxDepGenOrder','default',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},...
    {@()isa(params.R0,'polyZonotope')});
add2options('polyZono.maxPolyZonoRatio','default',{@isscalar,@isnumeric,@(val)val>0},...
    {'isscalar','isnumeric','gezero'},...
    {@()isa(params.R0,'polyZonotope')});
add2options('polyZono.restructureTechnique','default',{@ischar,...
    @(val)any(ismember(getMembers('restructureTechnique'),val))},...
    {'ischar','memberrestructureTechnique'},...
    {@()isa(params.R0,'polyZonotope')});

% lagrangeRem
add2options('lagrangeRem.simplify','default',{@ischar,...
    @(val)any(ismember(getMembers('lagrangeRem.simplify'),val))},...
    {'ischar','memberlagrangeRem.simplify'});
add2options('lagrangeRem.method','default',{@ischar,...
    @(val)any(ismember(getMembers('lagrangeRem.method'),val))},...
    {'ischar','memberlagrangeRem.method'});
add2options('lagrangeRem.tensorParallel','default',{@isscalar,@islogical},{'isscalar','islogical'});
add2options('lagrangeRem.replacements','optional',{@(val)isa(val,'function_handle')},{'isafunction_handle'});
add2options('lagrangeRem.zooMethods','mandatory',...
    {@ischar,@(val)any(ismember(getMembers('lagrangeRem.zooMethods'),val))},...
    {'ischar','memberlagrangeRem.zooMethods'},...
    {@()strcmp(options.lagrangeRem.method,'zoo')});
add2options('lagrangeRem.optMethod','default',...
    {@ischar,@(val)any(ismember(getMembers('lagrangeRem.optMethod'),val))},...
    {'ischar','memberlagrangeRem.optMethod'},...
    {@()strcmp(options.lagrangeRem.method,'taylorModel')});
add2options('lagrangeRem.maxOrder','optional',...
    {@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'},...
    {@()any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}))});
add2options('lagrangeRem.tolerance','optional',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},...
    {@()any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}))});
add2options('lagrangeRem.eps','optional',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},...
    {@()any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}))});

% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------

