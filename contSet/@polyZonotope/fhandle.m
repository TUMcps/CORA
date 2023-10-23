function [f,fvec,fGI] = fhandle(pZ,varargin)
% fhandle - returns the function handle for a given polyZonotope
%
% Syntax:
%    [f,fvec,fGI] = fhandle(pZ)
%    [f,fvec,fGI] = fhandle(pZ,ids)
%
% Inputs:
%    pZ - polyZonotope object
%    ids - cell array of ids
%
% Outputs:
%    f - function handle for dependent generators (+ center), where f has
%             length(ids) number of arguments
%    fvec - same as f, but only has one argument sorted according to
%             vertcat(ids{:})
%    fGI - function handle for rest zonotope
%
% Example: 
%    pZ = polyZonotope([0;0;0],diag([1,2,3]),zeros(3,0),eye(3),[1;2;3]);
%    f = fhandle(pZ,{[3;1],2});
%    f([1;0],2) % ans = [0;4;3];   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: subs

% Authors:       Victor Gassmann
% Written:       12-January-2021
% Last update:   11-March-2022
%                04-July-2022 (VG, update input check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
elseif nargin == 1
    id = pZ.id;
elseif nargin == 2
    if iscell(varargin{1})
%         for i=1:length(varargin{1})
%             inputArgsCheck({{varargin{1}{i},'att','double',{'integer','ncols',1}}});
%         end
        id = vertcat(varargin{1}{:});
    else
%         inputArgsCheck({{varargin{1},'att','double',{'integer','ncols',1}}});
        id = varargin{1};
    end
end

% make sure all ids are unique
if length(unique(id))<length(id)
    throw(CORAerror('CORA:wrongValue','second','ids need to be unique!'));
end
% make sure id contains all ids of pZ (and only those)
if ~isempty(setdiff(id,pZ.id)) || ~isempty(setdiff(pZ.id,id))
    throw(CORAerror('CORA:wrongValue','second','ids need to contain all ids!'));
end

% extract center
fc = pZ.c;
% bring to some order (by sorting)
[~,indpz] = sort(pZ.id);
[~,indid] = sort(id);

if isempty(pZ.G)
    fGt = @(x) fc;
else
    % compute function handle
    fGt = @(x) fc + sum(pZ.G.*prod(x(indid).^pZ.E(indpz,:),1),2);
end

% return from single id vector to number of specified inputs for function
% handle (as many as in original cell array)
fG = @(varargin) fGt(vertcat(varargin{:}));
f = @(varargin) fG(varargin{:}); 
% return alternatively single-id-vector version
fvec = @(x) fGt(x);
% return function handle for rest zonotope
fGI = @(r) pZ.GI*r;

% ------------------------------ END OF CODE ------------------------------
