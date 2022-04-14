function [paramsList,optionsList] = config_contDynamics_simulateNormal(sys,params,options)
% config_contDynamics_simulateNormal - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_contDynamics_simulateNormal(sys,params,options)
%
% Inputs:
%    sys - contDynamics object
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
add2params('U','default',{@(val)any(ismember(getMembers('Usim'),class(val)))},{'memberUsim'});
add2params('u','default',{@isnumeric},{'isnumeric'});
add2params('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
add2params('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});
add2params('paramInt','mandatory',{@(val)length(val)==sys.nrOfParam,...
    @(val)isa(val,'interval') || (isvector(val) && isnumeric(val))},...
    {'eqparam','vectororinterval'},{@()isa(sys,'nonlinParamSys')});
add2params('y0guess','mandatory',{@(val)length(val)==sys.nrOfConstraints},...
    {'eqconstr'},{@()isa(sys,'nonlinDASys')});

% append entries to list of algorithm parameters
add2options('points','mandatory',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,1)},...
    {'isscalar','isnumeric','integer','geone'});
add2options('fracVert','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'});
add2options('fracInpVert','mandatory',{@isscalar,@isnumeric,@(val)ge(val,0) && le(val,1)},...
    {'isscalar','isnumeric','normalized'});
add2options('inpChanges','mandatory',{@isscalar,@isnumeric,@(val)mod(val,1)==0,@(val)ge(val,0)},...
    {'isscalar','isnumeric','integer','gezero'});

% 3. prepare lists for output args
[paramsList,optionsList] = outputParamsOptionsLists();

end

%------------- END OF CODE --------------
