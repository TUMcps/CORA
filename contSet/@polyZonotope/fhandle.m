function [f,fvec,fGrest] = fhandle(pZ,ids)
% fhandle - returns the function handle for a given polyZonotope
%
% Syntax:  
%    [f,fvec,fGrest] = fhandle(pZ)
%    [f,fvec,fGrest] = fhandle(pZ,ids)
%
% Inputs:
%    pZ     - polynomial zonotope 
%    ids    - cell array of ids
%
% Outputs:
%    f      - function handle for dependent generators (+ center), where f has
%             length(ids) number of arguments
%    fvec   - same as f, but only has one argument sorted according to
%             vertcat(ids{:})
%    fGrest - function handle for rest zonotope
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

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if nargin==1
    id = pZ.id;
elseif nargin==2
    id = vertcat(ids{:});
    if numel(unique(id))~=numel(id)
        error('No ids can occur twice');
    end
    if numel(pZ.id)~=numel(id) || ~isempty(setdiff(pZ.id,id))
        error('pZ.id and ids must contain the same ids');
    end
end
fc = pZ.c;
[~,indpz] = sort(pZ.id);
[~,indid] = sort(id);

if isempty(pZ.G)
    fGt = @(x) fc;
else
    fGt = @(x) fc + sum(pZ.G.*prod(x(indid).^pZ.expMat(indpz,:),1),2);
end
fG = @(varargin) fGt(vertcat(varargin{:}));
f = @(varargin) fG(varargin{:}); 
fvec = @(x) fGt(x);
fGrest = @(r) pZ.Grest*r;