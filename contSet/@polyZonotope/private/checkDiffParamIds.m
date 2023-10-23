function [id_diff,id_param] = checkDiffParamIds(pZ,varargin)
% checkDiffParamIds - ???
%
% Syntax:
%    [id_diff,id_param] = checkDiffParamIds(pZ,varargin)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    id_diff - ???
%    id_param - ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(varargin)
    id_diff = pZ.id;
    id_param = zeros(0,1);
elseif length(varargin)==2
    id_diff = varargin{1};
    id_param = varargin{2};
    if all(size(id_param)==0)
        id_param = zeros(0,1);
    end
    % check input args
    inputArgsCheck({{id_diff,'att','double',{'ncols',1,'integer'}};
                    {id_param,'att','double',{'ncols',1,'integer'}}});
else
    throw(CORAerror('CORA:notSupported',...
        'Either one or three arguments have to be supplied!'));
end

% check if all ids are unique
if length(unique([id_diff;id_param]))<length([id_diff;id_param])
    throw(CORAerror('CORA:wrongValue','second/third',...
        'Ids need to be unique!'));
end

% check if id_diff + id_param makes up all of the ids
if ~isempty(setdiff(pZ.id,[id_diff;id_param]))
    throw(CORAerror('CORA:wrongValue','second/third',...
        'Ids need to contain all ids of pZ!'));
end

% ------------------------------ END OF CODE ------------------------------
