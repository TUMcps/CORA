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
% Last revision:---

%------------- BEGIN CODE --------------

% 1. init lists
initParamsOptionsLists();

% append entries to list of model parameters
add2params('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
add2params('U','default',{@(val)any(ismember(getMembers('U'),class(val)))},...
    {'memberU'},{@()isa(sys,'nonlinearSysDT')});
add2params('tStart','default',{@isscalar,@(val)ge(val,0)},...
    {'isscalar','gezero'});
add2params('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},...
    {'isscalar','getStart'});
add2params('u','default',{@isnumeric,@(val)c_inputTrajDT(val,sys,params)},{'isnumeric',''});
add2params('y','default',{@isnumeric,@(val)c_measurements(val,sys,params)},{'isnumeric',''});
% add2params('V','mandatory',{@(val)any(ismember(getMembers('V'),class(val))),...
%     @(val)eq(dim(val),sys.nrOfOutputs)},{'memberV','eqoutput'});
add2params('V','mandatory',{@(val)any(ismember(getMembers('V'),class(val)))},{'memberV'});
add2params('W','mandatory',{@(val)any(ismember(getMembers('W'),class(val))),...
    @(val)eq(dim(val),sys.dim)},{'memberW','eqsysdim'});


% append entries to list of algorithm parameters
add2options('verbose','default',{@isscalar,@islogical},...
    {'isscalar','islogical'});
add2options('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});
add2options('saveOrder','optional',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});

add2options('linAlg','default',{@ischar,@(val)any(ismember(getMembers('linAlg'),val))},...
    {'ischar','memberlinAlg'});
add2options('alg','mandatory',{@ischar,@(val)any(ismember(getMembers('alg4observe'),val))},...
    {'ischar','memberalg4observe'});
add2options('zonotopeOrder','mandatory',{@isscalar,@(val)ge(val,1)},...
    {'isscalar','geone'});

add2options('timeStep','mandatory',{@isscalar,@(val)val>0,...
    @(val)abs(params.tFinal/val - round(params.tFinal/val))<1e-9},...
    {'isscalar','gezero','intsteps'});

% for nonlinearSysDT
add2options('tensorOrder','mandatory',...
    {@isscalar,@(val)mod(val,1)==0,@(val)any(val==[2,3])},{'isscalar','integer','2or3'},{@()isa(sys,'nonlinearSysDT')});
add2options('errorOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'},{@()isa(sys,'nonlinearSysDT')});

% for simulateGaussian call
add2options('points','optional',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});
add2options('p_conf','optional',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'});

% for special solvers
add2options('solver','optional',{@ischar},{'ischar'});


% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------
