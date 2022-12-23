function [res,vars] = pre_volume(classname,S,varargin)
% pre_volume - pre-processing for contSet/volume functions including
%    1. setting of default values
%    2. checking of input arguments
%    3. computation of special cases (e.g., empty set)
%    the class name of the calling function is handed as an input argument
%    as a read-out from the call stack would be computationally inefficient
%
% Syntax:
%    [res,vars] = pre_volume(classname,S)
%    [res,vars] = pre_volume(classname,S,method)
%    [res,vars] = pre_volume(classname,S,method,order)
%
% Inputs:
%    classname - class of calling function
%    S - contSet object
%    method - (optional) method for evaluation
%    order - (optional) zonotope order
%
% Outputs:
%    res - true/false whether result has been computed here
%    vars - cell-array whose content depends on value of res:
%           res = true: value of volume
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
if nargin < 1
    throw(CORAerror('CORA:notEnoughInputArgs',1));
elseif nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs',4));
end

% check input arguments
inputArgsCheck({{S,'att',classname}});

% check input arguments: two tries depending on value for type
if strcmp(classname,'zonotope')
    % only zonotopes: set default values method ('exact') and order (5)
    [method,order] = setDefaultValues({'exact',5},varargin);
    % check input arguments
    inputArgsCheck({{method,'str',{'exact','reduce','alamo'}},...
                    {order,'att','numeric',{'integer','positive'}}});
end

% special case: contSet is empty set
if isemptyobject(S)
    res = true;
    vars{1} = 0; return
end

% result has to be computed in calling function
res = false;
vars{1} = S;
if strcmp(classname,'zonotope')
    vars{2} = method;
    vars{3} = order;
end

%------------- END OF CODE --------------
