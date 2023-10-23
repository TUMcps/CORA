function print(pZ,varargin)
% print - prints a polynomial zonotope on the command window
%
% Syntax:
%    print(pZ)
%    print(pZ,'Digits',digits)
%    print(pZ,'Ids',ids,'Vars',vars)
%    print(pZ,'Ids',ids,'Vars',vars,'Digits',digits)
%
% Inputs:
%    pZ - polyZonotope object
%    Name-Value pairs (all options, arbitrary order):
%       <'Digits',digits> - number of printed digits
%       <'Ids',ids> - cell array where each vector member represents ids to
%                     be printed with char given in chars
%       <'Vars',vars> - cell array of same length as ids containing
%                         variable "names"
%
% Outputs:
%    -
%
% Example: 
%    pZ = polyZonotope([pi;-1],[2,0;1,1],zeros(2,0),[1,2;2,3;0,1],[1;2;3]);
%    print(pZ,'Ids',{[1;3],2},'Vars',{'x','y'},'Digits',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       12-January-2021 
% Last update:   07-June-2022 (MW, switch to name-value syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always uneven
if mod(nargin,2) ~= 1
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check for redundant name-value pairs
    checkNameValuePairs(NVpairs,{'Digits','Ids','Vars'});
    % dimension given?
    [NVpairs,digits] = readNameValuePair(NVpairs,'Digits','isscalar',2);
    % center given?
    [NVpairs,ids] = readNameValuePair(NVpairs,'Ids',{},{pZ.id});
    % number of generators given?
    [NVpairs,vars] = readNameValuePair(NVpairs,'Vars',{},{'b'});
end

% compute id
id = vertcat(ids{:});

% check input arguments
if length(ids)~=length(vars)
     throw(CORAerror('CORA:wrongValue',"name-value pair 'Ids','Vars'",...
         "Values for 'Ids' and 'Vars' have to have same length."));
end
if numel(pZ.id) ~= numel(id) || ~isempty(setdiff(pZ.id,id))
    throw(CORAerror('CORA:wrongValue',"name-value pair 'Ids'",...
        "Vector for 'Ids' should contain same identifiers as in polynomial zonotope."));
end


[~,f] = fhandle(pZ,ids);
x = sym('x',[length(id),1]);
% fill symbolic placeholder with variables of specified name
for i=1:length(ids)
    k = numel(vertcat(ids{1:(i-1)}));
    x(k+1:k+numel(ids{i})) = sym(vars{i},size(ids{i}));
end
res = f(x);

% check if pZ has GI
if ~isempty(pZ.GI) && ~all(pZ.GI==0,'all')
    ind_zero = all(pZ.GI==0,2);
    Zrest = sym('Zrest',[sum(~ind_zero),1]);
    res(~ind_zero) = res(~ind_zero) + Zrest;
end

% display object
disp(vpa(res,digits));

% ------------------------------ END OF CODE ------------------------------
