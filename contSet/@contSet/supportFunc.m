function varargout = supportFunc(S,dir,varargin)
% supportFunc - evaluates the support function of a set along a given direction
%
% Description:
%    computes \max_{x \in \mathcal{S}}  l^T * x
%
% Syntax:
%    [val,x,fac] = supportFunc(S)
%    [val,x,fac] = supportFunc(S,dir)
%    [val,x,fac] = supportFunc(S,dir,type)
%    [val,x,fac] = supportFunc(S,dir,type,method)
%    [val,x,fac] = supportFunc(S,dir,type,method,vals)
%    [val,x,fac] = supportFunc(S,dir,type,method,vals,tol)
%
% Inputs:
%    S - contSet object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - minimum ('lower'), maximum ('upper') or range ('range')
%    method - only polyZonotope/conPolyZonotope: method that is used to
%             calculate the bounds for the dependent part of the polynomial
%             zonotope
%    splits - number of splits that are performed to calculate the bounds
%    maxOrder - maximum polynomial order of the Taylor model
%    tol - tolerance
%
% Outputs:
%    val - value of support function
%    x - support vector
%    fac - factors for set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       18-August-2022
% Last update:   23-November-2022 (MW, add classname as input argument)
% Last revision: 27-March-2023 (MW, restructure relation to subclass)

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs',2));
elseif nargin > 6
    throw(CORAerror('CORA:tooManyInputArgs',6));
end

% parse input arguments
[type,method,maxOrderOrSplits,tol] = ...
    setDefaultValues({'upper','interval',8,1e-3},varargin);

if isa(S,'polyZonotope')
    % polyZonotope: set, direction, type, method, vals, tol    

    % check input arguments
    inputArgsCheck({{S,'att','contSet','scalar'};
                    {dir,'att','numeric','vector'};
                    {type,'str',{'lower','upper','range'}};
                    {method,'str',{'interval','split','bnb','bnbAdv',...
                        'globOpt','bernstein','quadProg'}}; ...
                    {maxOrderOrSplits,'att','numeric',{'scalar','integer','positive'}};...
                    {tol,'att','numeric',{'nonnegative','scalar'}}});

    % meaning of maxOrderOrSplits depends on method

elseif isa(S,'conPolyZono')

    % check input arguments
    inputArgsCheck({{S,'att','contSet','scalar'};
                    {dir,'att','numeric','vector'};
                    {type,'str',{'lower','upper','range'}};
                    {method,'str',{'interval','split','conZonotope','quadProg'}};
                    {maxOrderOrSplits,'att','numeric',{'scalar','integer','positive'}}});

else

    % check input arguments
    inputArgsCheck({{S,'att','contSet'};
                    {dir,'att','numeric','vector'};
                    {type,'str',{'upper','lower','range'}}});
    
end

% transpose dir if necessary
if size(dir,1) == 1 && size(dir,2) > 1
    dir = dir';
end

% ensure that dimension of direction fits dimension of the set
if dim(S) ~= length(dir)
    throw(CORAerror('CORA:wrongValue','second',...
        [num2str(dim(S)) '-dimensional column vector.']));
end

% result has to be computed in calling function
varargout = cell(1,max(nargout,1));
try
    [varargout{:}] = supportFunc_(S,dir,type,method,maxOrderOrSplits,tol);
catch ME
    % empty set case
    if representsa_(S,'emptySet',eps)
        if strcmp(type,'upper')
            varargout{1} = -Inf;
        elseif strcmp(type,'lower')
            varargout{1} = +Inf;
        elseif strcmp(type,'range')
            varargout{1} = [-Inf,+Inf];
        end
        varargout{2} = [];
        varargout{3} = [];
    else
        rethrow(ME);
    end

end

% ------------------------------ END OF CODE ------------------------------
