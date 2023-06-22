function [paramsList,optionsList] = config_contDynamics_observe(sys,params,options)
% config_contDynamics_observe - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_contDynamics_observe(sys,params,options)
%
% Inputs:
%    sys - linearSysDT or nonlinearSysDT object
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

% Author:       Mark Wetzlinger, Matthias Althoff
% Written:      06-July-2021
% Last update:  07-July-2021 (MA, input set U removed)
% Last revision:19-June-2023 (MW, structs, remove global variables)

%------------- BEGIN CODE --------------

% list of model parameters
paramsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of model parameters
paramsList(end+1,1) = add2list('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
paramsList(end+1,1) = add2list('U','default',{@(val)any(ismember(getMembers('U'),class(val)))},...
    {'memberU'},{@()isa(sys,'nonlinearSysDT')});
paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@(val)ge(val,0)},...
    {'isscalar','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},...
    {'isscalar','getStart'});
paramsList(end+1,1) = add2list('u','default',{@isnumeric,@(val)c_inputTrajDT(val,sys,params)},{'isnumeric',''});
paramsList(end+1,1) = add2list('y','default',{@isnumeric,@(val)c_measurements(val,sys,params)},{'isnumeric',''});
paramsList(end+1,1) = add2list('V','mandatory',{@(val)any(ismember(getMembers('V'),class(val))),...
    @(val)c_V(val,sys,params,options)},{'memberV',''});
paramsList(end+1,1) = add2list('W','mandatory',{@(val)any(ismember(getMembers('W'),class(val))),...
    @(val)eq(dim(val),sys.dim)},{'memberW','eqsysdim'});


% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of algorithm parameters
optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},...
    {'isscalar','islogical'});
optionsList(end+1,1) = add2list('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
optionsList(end+1,1) = add2list('saveOrder','optional',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});

optionsList(end+1,1) = add2list('linAlg','default',{@ischar,@(val)any(ismember(getMembers('linAlg'),val))},...
    {'ischar','memberlinAlg'});
optionsList(end+1,1) = add2list('alg','mandatory',{@ischar,@(val)any(ismember(getMembers('alg4observe'),val))},...
    {'ischar','memberalg4observe'});
optionsList(end+1,1) = add2list('zonotopeOrder','mandatory',{@isscalar,@(val)ge(val,1)},...
    {'isscalar','geone'});

optionsList(end+1,1) = add2list('timeStep','mandatory',{@isscalar,@(val)val>0,...
    @(val)abs(params.tFinal/val - round(params.tFinal/val))<1e-9},...
    {'isscalar','gezero','intsteps'});

% for nonlinearSysDT
optionsList(end+1,1) = add2list('tensorOrder','mandatory',...
    {@isscalar,@(val)mod(val,1)==0,@(val)any(val==[2,3])},{'isscalar','integer','2or3'},{@()isa(sys,'nonlinearSysDT')});
optionsList(end+1,1) = add2list('errorOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'},{@()isa(sys,'nonlinearSysDT')});

% for simulateGaussian call
optionsList(end+1,1) = add2list('points','optional',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});
optionsList(end+1,1) = add2list('p_conf','optional',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'});

% for special solvers
optionsList(end+1,1) = add2list('solver','optional',{@ischar},{'ischar'});

end

%------------- END OF CODE --------------
