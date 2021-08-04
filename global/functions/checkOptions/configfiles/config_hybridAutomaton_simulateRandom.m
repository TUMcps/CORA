function [paramsList,optionsList] = config_hybridAutomaton_simulateRandom(sys,params,options)
% config_hybridAutomaton_simulateRandom - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_hybridAutomaton_simulateRandom(sys,params,options)
%
% Inputs:
%    sys - hybridAutomaton object
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
% Written:      04-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. init lists
initParamsOptionsLists();

% append entries to list of model parameters
add2params('tStart','default',{@isscalar,@isnumeric,@(val)ge(val,0)},{'isscalar','isnumeric','gezero'});
add2params('tFinal','mandatory',{@isscalar,@isnumeric,@(val)val>0},{'isscalar','isnumeric','gezero'});
add2params('startLoc','mandatory',{@isscalar,@isnumeric,...
    @(val)mod(val,1)==0,@(val)ge(val,1),@(val)le(val,length(sys.location))},...
    {'isscalar','isnumeric','integer','geone','leloc'});
add2params('finalLoc','default',{@isscalar,@isnumeric,...
    @(val)mod(val,1)==0,@(val)ge(val,0),@(val)le(val,length(sys.location))},...
    {'isscalar','isnumeric','integer','geone','leloc'});
add2params('R0','mandatory',{@(val)isa(val,'contSet'),...
    @(val)dim(val)==sys.location{params.startLoc}.contDynamics.dim},{'isacontSet','eqsysdim'});
add2params('U','default',{@(val) (~iscell(val) && isa(val,'contSet') ) || ...
    ( iscell(val) && (all(size(val) == [length(sys.location),1]) || ...
        all(size(val) == [1,length(sys.location)])) )},{'isacontSet','???'});
    
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


