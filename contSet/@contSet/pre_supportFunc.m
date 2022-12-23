function [res,vars] = pre_supportFunc(classname,S,dir,varargin)
% pre_supportFunc - pre-processing for contSet/supportFunc functions including
%    1. setting of default values
%    2. checking of input arguments
%    3. computation of special cases (e.g., empty set)
%    the class name of the calling function is handed as an input argument
%    as a read-out from the call stack would be computationally inefficient
%
% Syntax:
%    [res,vars] = pre_supportFunc(classname,S)
%    [res,vars] = pre_supportFunc(classname,S,dir)
%    [res,vars] = pre_supportFunc(classname,S,dir,type)
%    [res,vars] = pre_supportFunc(classname,S,dir,type,method)
%    [res,vars] = pre_supportFunc(classname,S,dir,type,method,vals)
%    [res,vars] = pre_supportFunc(classname,S,dir,type,method,vals,tol)
%
% Inputs:
%    classname - class of calling function
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
%    res - true/false whether result has been computed here
%    vars - cell-array whose content depends on value of res:
%           res = true: value of supportFunc
%           res = false: checked input argument arguments incl. default
%                        values for volume function evaluation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      18-August-2022
% Last update:  23-November-2022 (MW, add classname as input argument)
% Last revision:---

%------------- BEGIN CODE --------------

% check number of input arguments
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs',2));
elseif nargin > 7
    throw(CORAerror('CORA:tooManyInputArgs',7));
end

% parse input arguments
[type,method,maxOrderOrSplits,tol] = ...
    setDefaultValues({'upper','interval',8,1e-3},varargin);

if strcmp(classname,'polyZonotope')
    % polyZonotope: set, direction, type, method, vals, tol    

    % check input arguments
    inputArgsCheck({{S,'att',classname,'scalar'};
                    {dir,'att','numeric','vector'};
                    {type,'str',{'lower','upper','range'}};
                    {method,'str',{'interval','split','bnb','bnbAdv',...
                        'globOpt','bernstein','quadProg'}}; ...
                    {maxOrderOrSplits,'att','numeric',{'scalar','integer','positive'}};...
                    {tol,'att','numeric',{'nonnegative','scalar'}}});

    % meaning of maxOrderOrSplits depends on method

elseif strcmp(classname,'conPolyZono')

    % check input arguments
    inputArgsCheck({{S,'att',classname,'scalar'};
                    {dir,'att','numeric','vector'};
                    {type,'str',{'lower','upper','range'}};
                    {method,'str',{'interval','split','conZonotope','quadProg'}};
                    {maxOrderOrSplits,'att','numeric',{'scalar','integer','positive'}}});

else

    % check input arguments
    inputArgsCheck({{S,'att',classname};
                    {dir,'att','numeric','vector'};
                    {type,'str',{'upper','lower','range'}}});
    
end


% empty set case
if isemptyobject(S)
    res = true;
    if strcmp(type,'upper')
        vars{1} = -Inf;
    elseif strcmp(type,'lower')
        vars{1} = +Inf;
    elseif strcmp(type,'range')
        vars{1} = [-Inf,+Inf];
    end
    return
end

% transpose dir if necessary
if size(dir,1) == 1 && size(dir,2) >= 1
    dir = dir';
end

% result has to be computed in calling function
res = false;
vars = {S,dir,type,method,maxOrderOrSplits,tol};

%------------- END OF CODE --------------
