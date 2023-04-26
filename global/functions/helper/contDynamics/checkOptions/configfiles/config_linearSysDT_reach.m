function [paramsList,optionsList] = config_linearSysDT_reach(sys,params,options)
% config_linearSysDT_reach - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_linearSysDT_reach(sys,params,options)
%
% Inputs:
%    sys - hybridSystem or contDynamics object
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
% Written:      25-January-2021
% Last update:  03-February-2021
% Last revision:---

%------------- BEGIN CODE --------------

% 1. init lists
initParamsOptionsLists();

% append entries to list of model parameters
add2params('tStart','default',{@isscalar,@(val)ge(val,0)},...
    {'isscalar','gezero'});
add2params('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},...
    {'isscalar','getStart'});
add2params('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),...
    @(val)eq(dim(val),sys.dim)},{'memberR0','eqsysdim'});
add2params('U','default',{@(val)any(ismember(getMembers('U'),class(val))),...
    @(val)eq(dim(val),sys.nrOfInputs)},{'memberU','eqinput'});
add2params('u','default',{@isnumeric,@(val)eq(size(val,1),sys.nrOfInputs),...
    @(val)c_inputTrajDT(val,sys,params)},{'isnumeric','eqinput',''});
add2params('W','default',{@(val)any(ismember(getMembers('W'),class(val))),...
    @(val)eq(dim(val),sys.dim)},{'memberW','eqsysdim'});
add2params('V','default',{@(val)any(ismember(getMembers('V'),class(val))),...
    @(val)eq(dim(val),sys.nrOfOutputs)},{'memberV','eqoutput'});

% append entries to list of algorithm parameters
add2options('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});
add2options('saveOrder','optional',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});
add2options('zonotopeOrder','mandatory',{@isscalar,@(val)ge(val,1)},{'isscalar','geone'});
add2options('reductionTechnique','default',...
    {@ischar,@(val)any(ismember(getMembers('reductionTechnique'),val))},...
    {'ischar','memberreductionTechnique'});

add2options('compOutputSet','default',{@isscalar,@islogical},...
    {'isscalar','islogical'});

% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------

