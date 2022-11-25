function print(pZ,varargin)
% print - print polynomial zonotope
%
% Syntax:  
%    print(pZ)
%    print(pZ,digits)
%    print(pZ,ids,chars)
%    print(pZ,ids,chars,digits)
%
% Inputs:
%    pZ     - polyZonotope to print
%    digits - number of printed digits
%    ids    - cell array where each vector member represents ids to be
%             printed with char given in chars
%    chars  - cell array of same length as ids containing variable "names"
%
% Example: 
%    pZ = polyZonotope([pi;-1],[2,0;1,1],zeros(2,0),[1,2;2,3;0,1],[1;2;3]);
%    print(pZ,{[1;3],2},{'x','y'},3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Victor Gassmann
% Written:      12-January-2021 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------
digits = 2;
if nargin==0
    error('Argument required');
end
if nargin<=2
    ids = {pZ.id};
    vars = {'b'};
    id = pZ.id;
    if nargin==2
        digits = varargin{1};
    end
elseif nargin<=4
    ids = varargin{1};
    vars = varargin{2};
    if length(ids)~=length(vars)
        error('ids and vars have to have same length');
    end
    id = vertcat(ids{:});
    if numel(pZ.id)~=numel(id) || ~isempty(setdiff(pZ.id,id))
        error('pZ.id and ids must contain the same ids');
    end
    if nargin==4
        digits = varargin{3};
    end
else
    error('Wrong number of arguments');
end
[~,f] = fhandle(pZ,ids);
x = sym('x',[length(id),1]);
% fill symbolic placeholder with variables of specified name
for i=1:length(ids)
    k = numel(vertcat(ids{1:(i-1)}));
    x(k+1:k+numel(ids{i})) = sym(vars{i},size(ids{i}));
end
res = f(x);
% check if pZ has Grest
if ~isempty(pZ.Grest) && ~all(pZ.Grest==0,'all')
    ind_zero = all(pZ.Grest==0,2);
    Zrest = sym('Zrest',[sum(~ind_zero),1]);
    res(~ind_zero) = res(~ind_zero) + Zrest;
end
disp(vpa(res,digits));
%------------- END OF CODE --------------