function [paramsList,optionsList] = config_linProbSys_reach(sys,params,options)
% config_linProbSys_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_linProbSys_reach(sys,params,options)
%
% Inputs:
%    sys - linParamSys object
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
paramsList(end+1,1) = add2list('u','default',{@isnumeric},{'isnumeric'});
paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of algorithm parameters
optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('compTimePoint','default',{@isscalar,@islogical},{'isscalar','islogical'});
optionsList(end+1,1) = add2list('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
optionsList(end+1,1) = add2list('saveOrder','optional',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});
optionsList(end+1,1) = add2list('gamma','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0)},{'isscalar','isnumeric','gezero'});

optionsList(end+1,1) = add2list('timeStep','mandatory',{@isscalar,@(val)val>0,...
    @(val)abs(params.tFinal/val - round(params.tFinal/val))<1e-9,...
	@(val)c_inputTraj(val,sys,params,options)},...
    {'isscalar','gezero','intsteps',''});
optionsList(end+1,1) = add2list('taylorTerms','mandatory',{@isscalar,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','integer','geone'});
optionsList(end+1,1) = add2list('zonotopeOrder','mandatory',{@isscalar,@isnumeric,@(val)ge(val,1)},...
    {'isscalar','isnumeric','geone'});

optionsList(end+1,1) = add2list('compOutputSet','default',{@isscalar,@islogical},...
    {'isscalar','islogical'});

end

%------------- END OF CODE --------------
