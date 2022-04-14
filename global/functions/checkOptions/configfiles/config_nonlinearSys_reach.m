function [paramsList,optionsList] = config_nonlinearSys_reach(sys,params,options)
% config_nonlinearSys_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_nonlinearSys_reach(sys,params,options)
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
% Written:      26-January-2021
% Last update:  02-February-2021
% Last revision:---

%------------- BEGIN CODE --------------

% 1. init lists
initParamsOptionsLists();

% append entries to list of model parameters
add2params('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
add2params('U','default',{@(val)any(ismember(getMembers('U'),class(val)))},{'memberU'});
add2params('u','default',{@isnumeric},{'isnumeric'});
add2params('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
add2params('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});

% append entries to list of algorithm parameters
add2options('timeStep','mandatory',{@isscalar,@(val)val>0,...
    @(val)abs(params.tFinal/val - round(params.tFinal/val))<1e-9,...
	@(val)c_inputTraj(val,sys,params,options)},...
    {'isscalar','gezero','intsteps',''},{@()~contains(options.alg,'adaptive')});
add2options('taylorTerms','mandatory',{@isscalar,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','integer','geone'},{@()~contains(options.alg,'adaptive')});

add2options('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
add2options('saveOrder','optional',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});

add2options('zonotopeOrder','mandatory',{@isscalar,@(val)ge(val,1)},...
    {'isscalar','geone'},{@()~contains(options.alg,'adaptive')});
add2options('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique4nlsys'),val))},...
    {'ischar','memberreductionTechnique4nlsys'});

% only dep. factor reachSet options (parameters/options used from AROC)
add2options('approxDepOnly','optional',{@isscalar,@islogical},{'isscalar','islogical'});
add2options('POpt.j','optional',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)val>=1},...
    {'isscalar','isnumeric','integer','geone'});
add2options('POpt.idv','optional',{@(val)true},{'none'});
add2options('POpt.idp','optional',{@(val)true},{'none'});
add2options('POpt.np','optional',{@(val)true},{'none'});
add2options('POpt.Q','optional',{@(val)true},{'none'});
add2options('POpt.reachOpt','optional',{@(val)true},{'none'});
add2options('POpt.pc_1j','optional',{@(val)true},{'none'});
add2options('POpt.D_1j','optional',{@(val)true},{'none'});
add2options('POpt.idp_1j','optional',{@(val)true},{'none'});
add2options('POpt.A_pz_1j','optional',{@(val)true},{'none'});
add2options('POpt.b_pz_1j','optional',{@(val)true},{'none'});
add2params('refPoints','optional',{@isnumeric,@(val)size(val,1)==sys.dim,...
    @(val)size(val,2)==reachSteps(params,options)+1},...
    {'isnumeric','eqsysdim','eqreachSteps'});

add2options('reductionInterval','default',{@isscalar,@(val)ge(val,1),@(val)isinf(val)||mod(val,1)==0},...
    {'isscalar','geone','integerorInf'});
add2options('maxError','default',{@isvector,@(val)all(ge(val,0)),@(val)length(val)==sys.dim},...
    {'isvector','vectorgezero','eqsysdim'});

add2options('alg','mandatory',{@ischar,@(val)any(ismember(getMembers('alg'),val))},...
    {'ischar','memberalg'});
add2options('tensorOrder','mandatory',...
    {@isscalar,@(val)mod(val,1)==0,@(val)c_tensorOrder(val,sys,options)},...
    {'isscalar','integer',''},{@()~contains(options.alg,'adaptive')});

add2options('errorOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'},...
    {@()~contains(options.alg,'adaptive') && options.tensorOrder>2});
add2options('errorOrder3','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'},...
    {@()~contains(options.alg,'adaptive') && options.tensorOrder>3});
add2options('intermediateOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'},...
    {@()~contains(options.alg,'adaptive') && ~(strcmp(options.alg,'lin') && options.tensorOrder == 2)});

% polyZono
add2options('polyZono.maxDepGenOrder','default',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'},...
    {@()contains(options.alg,'poly')});
add2options('polyZono.maxPolyZonoRatio','default',{@isscalar,@isnumeric,@(val)val>0},...
    {'isscalar','isnumeric','gezero'},...
    {@()contains(options.alg,'poly')});
add2options('polyZono.restructureTechnique','default',{@ischar,...
    @(val)any(ismember(getMembers('restructureTechnique'),val))},...
    {'ischar','memberrestructureTechnique'},...
    {@()contains(options.alg,'poly')});

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

