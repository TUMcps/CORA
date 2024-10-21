function [paramsList,optionsList] = config_linearSys_reachBackward
% config_linearSys_reachBackward - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_linearSys_reachBackward
%
% Inputs:
%    -
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

% Authors:       Mark Wetzlinger
% Written:       15-December-2022
% Last update:   12-July-2023 (MW, adapt to new syntax)
%                21-September-2024 (MW, adapt to new syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init structs
[paramsList,optionsList] = initDynParameterList();

% list of model parameters ------------------------------------------------

% mandatory
paramsList(end+1,1) = add2list('tFinal','mandatory');
paramsList(end+1,1) = add2list('R0','mandatory');

% default
paramsList(end+1,1) = add2list('tStart','default');
paramsList(end+1,1) = add2list('U','default');
paramsList(end+1,1) = add2list('u','default');
paramsList(end+1,1) = add2list('tu','default');
paramsList(end+1,1) = add2list('W','default');
paramsList(end+1,1) = add2list('V','default');


% optional

% list of algorithm parameters --------------------------------------------

% mandatory
optionsList(end+1,1) = add2list('timeStep','mandatory');
optionsList(end+1,1) = add2list('linAlg','mandatory');

% default
optionsList(end+1,1) = add2list('verbose','default');

% optional


%%% old lists (for reference)
% paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});
% paramsList(end+1,1) = add2list('Rend','mandatory',{@(val)any(ismember(getMembers('U'),class(val))),@(val)eq(dim(val),sys.dim)},...
%     {'memberRend','eqsysdim'});
% 
% paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
% paramsList(end+1,1) = add2list('U','default',{@(val)any(ismember(getMembers('U'),class(val))),...
%     @(val)eq(dim(val),sys.nrOfInputs)},{'memberU','eqinput'});
% paramsList(end+1,1) = add2list('u','default',{@isnumeric,@(val)eq(size(val,1),sys.nrOfInputs)},{'isnumeric','eqinput'});
% paramsList(end+1,1) = add2list('tu','default',{@isvector,@isnumeric,@(val)all(diff(val)>0),...
%     @(val)c_tu(val,sys,params,options)},{'isvector','isnumeric','vectorgezero',''},...
%     {@()isfield(options,'linAlg') && strcmp(options.linAlg,'adaptive')});
% paramsList(end+1,1) = add2list('W','default',{@(val)any(ismember(getMembers('W'),class(val))),...
%     @(val)eq(dim(val),sys.dim)},{'memberW','eqsysdim'});
% paramsList(end+1,1) = add2list('V','default',{@(val)any(ismember(getMembers('V'),class(val))),...
%     @(val)eq(dim(val),sys.nrOfOutputs)},{'memberV','eqoutput'});
% 
% 
% optionsList(end+1,1) = add2list('timeStep','mandatory',{@isscalar,@(val)val>0,...
%     @(val)withinTol((params.tFinal-params.tStart)/val,round((params.tFinal-params.tStart)/val),1e-9),...
% 	@(val)c_inputTraj(val,sys,params,options)},...
%     {'isscalar','gezero','intsteps',''});
% 
% optionsList(end+1,1) = add2list('linAlg','mandatory',{@ischar,...
%     @(val)any(ismember(getMembers('linAlg4backward'),val))},...
%     {'ischar','memberlinAlg4backward'});
% 
% optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});

% ------------------------------ END OF CODE ------------------------------
