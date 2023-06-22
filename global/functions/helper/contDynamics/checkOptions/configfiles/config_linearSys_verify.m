function [paramsList,optionsList] = config_linearSys_verify(sys,params,options)
% config_linearSys_verify - configuration file for validation of
%    model parameters and algorithm parameters
%
% Syntax:
%    [paramsList,optionsList] = config_linearSys_verify(sys,params,options)
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
% Written:      22-April-2022
% Last update:  ---
% Last revision:19-June-2023 (MW, structs, remove global variables)

%------------- BEGIN CODE --------------

% list of model parameters
paramsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of model parameters
paramsList(end+1,1) = add2list('tStart','default',{@isscalar,@(val)ge(val,0)},{'isscalar','gezero'});
paramsList(end+1,1) = add2list('tFinal','mandatory',{@isscalar,@(val)ge(val,params.tStart)},{'isscalar','getStart'});
paramsList(end+1,1) = add2list('R0','mandatory',{@(val)any(ismember(getMembers('R0'),class(val))),@(val)eq(dim(val),sys.dim)},...
    {'memberR0','eqsysdim'});
paramsList(end+1,1) = add2list('U','default',{@(val)any(ismember(getMembers('U'),class(val))),...
    @(val)eq(dim(val),sys.nrOfInputs)},{'memberU','eqinput'});
paramsList(end+1,1) = add2list('u','default',{@isnumeric,@(val)eq(size(val,1),sys.nrOfInputs)},{'isnumeric','eqinput'});
paramsList(end+1,1) = add2list('tu','default',{@isvector,@isnumeric,@(val)all(diff(val)>0),...
    @(val)c_tu(val,sys,params,options)},{'isvector','isnumeric','vectorgezero',''});
paramsList(end+1,1) = add2list('W','default',{@(val)any(ismember(getMembers('W'),class(val))),...
    @(val)eq(dim(val),sys.dim)},{'memberW','eqsysdim'});
paramsList(end+1,1) = add2list('V','default',{@(val)any(ismember(getMembers('V'),class(val))),...
    @(val)eq(dim(val),sys.nrOfOutputs)},{'memberV','eqoutput'});
paramsList(end+1,1) = add2list('safeSet','default',{@(val)c_safeSet(val,sys,params,options)},{''});
paramsList(end+1,1) = add2list('unsafeSet','default',{@(val)c_unsafeSet(val,sys,params,options)},{''});

% list of algorithm parameters
optionsList = struct('name',{},'status',{},'checkfun',{},'errmsg',{},'condfun',{});

% append entries to list of algorithm parameters
optionsList(end+1,1) = add2list('verbose','default',{@isscalar,@islogical},{'isscalar','islogical'});

end

%------------- END OF CODE --------------
