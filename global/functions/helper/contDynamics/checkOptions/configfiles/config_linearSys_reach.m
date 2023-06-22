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
% Last revision:19-June-2023 (MW, structs, remove global variables)

%------------- BEGIN CODE --------------

% list of model parameters
paramsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});
paramsList(end+1,1) = add2list('R0','mandatory',...
    {@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
paramsList(end+1,1) = add2list('U','default',{@(val)any(ismember(getMembers('U'),class(val))),...
    @(val)eq(dim(val),sys.nrOfInputs)},{'memberU','eqinput'});
paramsList(end+1,1) = add2list('u','default',{@isnumeric,@(val)eq(size(val,1),sys.nrOfInputs)},{'isnumeric','eqinput'});
paramsList(end+1,1) = add2list('tu','default',{@isvector,@isnumeric,@(val)all(diff(val)>0),...
    @(val)c_tu(val,sys,params,options)},{'isvector','isnumeric','vectorgezero',''},...
    {@()isfield(options,'linAlg') && strcmp(options.linAlg,'adaptive')});
paramsList(end+1,1) = add2list('W','default',{@(val)any(ismember(getMembers('W'),class(val))),...
    @(val)eq(dim(val),sys.dim)},{'memberW','eqsysdim'});
paramsList(end+1,1) = add2list('V','default',{@(val)any(ismember(getMembers('V'),class(val))),...
    @(val)eq(dim(val),sys.nrOfOutputs)},{'memberV','eqoutput'});

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
optionsList(end+1,1) = add2list('saveOrder','optional',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});

optionsList(end+1,1) = add2list('linAlg','default',{@ischar,@(val)any(ismember(getMembers('linAlg'),val))},...
    {'ischar','memberlinAlg'});
optionsList(end+1,1) = add2list('error','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},{@()strcmp(options.linAlg,'adaptive')});

optionsList(end+1,1) = add2list('timeStep','mandatory',{@isscalar,@(val)val>0,@(val)abs(params.tFinal/val - round(params.tFinal/val))<1e-9,...
	@(val)c_inputTraj(val,sys,params,options)},...
    {'isscalar','gezero','intsteps',''},{@()~strcmp(options.linAlg,'adaptive')});
optionsList(end+1,1) = add2list('taylorTerms','mandatory',{@isscalar,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','integer','geone'},{@()~strcmp(options.linAlg,'adaptive')});

optionsList(end+1,1) = add2list('zonotopeOrder','mandatory',{@isscalar,@(val)ge(val,1)},...
    {'isscalar','geone'},{@()~any(strcmp(options.linAlg,{'adaptive','supportFunc'}))});
% optionsList(end+1,1) = add2list('l','mandatory',{@isnumeric,@(val)eq(size(val,1),sys.nrOfOutputs)},...
%     {'isnumeric','eqoutput'},{@()strcmp(options.linAlg,'supportFunc')});

optionsList(end+1,1) = add2list('compOutputSet','default',{@isscalar,@islogical},...
    {'isscalar','islogical'});

optionsList(end+1,1) = add2list('partition','mandatory',{@(val)c_partition(val,sys,options)},...
    {''},{@()strcmp(options.linAlg,'decomp')});
optionsList(end+1,1) = add2list('krylovError','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0)},...
    {'isscalar','isnumeric','gezero'},{@()strcmp(options.linAlg,'krylov')});
optionsList(end+1,1) = add2list('krylovStep','mandatory',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'},{@()strcmp(options.linAlg,'krylov')});

end

%------------- END OF CODE --------------
