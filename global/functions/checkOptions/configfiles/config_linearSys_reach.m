function [paramsList,optionsList] = config_linearSys_reach(sys,params,options)
% config_linearSys_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_linearSys_reach(sys,params,options)
%
% Inputs:
%    sys - linearSys object
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
add2params('u','default',{@isnumeric},{'isnumeric'});
add2params('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
add2params('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});

% append entries to list of algorithm parameters
add2options('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
add2options('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
add2options('saveOrder','optional',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});

add2options('linAlg','default',{@ischar,@(val)any(ismember(getMembers('linAlg'),val))},...
    {'ischar','memberlinAlg'});
add2options('error','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},{@()strcmp(options.linAlg,'adaptive')});

add2options('timeStep','mandatory',{@isscalar,@(val)val>0,@(val)abs(params.tFinal/val - round(params.tFinal/val))<1e-9,...
	@(val)c_inputTraj(val,sys,params,options)},...
    {'isscalar','gezero','intsteps',''},{@()~strcmp(options.linAlg,'adaptive')});
add2options('taylorTerms','mandatory',{@isscalar,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','integer','geone'},{@()~strcmp(options.linAlg,'adaptive')});
add2options('zonotopeOrder','mandatory',{@isscalar,@(val)ge(val,1)},...
    {'isscalar','geone'},{@()~strcmp(options.linAlg,'adaptive')});

add2options('partition','mandatory',{@(val)c_partition(val,sys,options)},...
    {''},{@()strcmp(options.linAlg,'decomp')});
add2options('krylovError','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},{@()strcmp(options.linAlg,'krylov')});
add2options('krylovStep','mandatory',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'},{@()strcmp(options.linAlg,'krylov')});


% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------
