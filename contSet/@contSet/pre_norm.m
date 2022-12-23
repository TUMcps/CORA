function [res,vars] = pre_norm(classname,S,varargin)
% pre_norm - pre-processing for contSet/norm functions including
%    1. setting of default values
%    2. checking of input arguments
%    3. computation of special cases (e.g., empty set)
%    the class name of the calling function is handed as an input argument
%    as a read-out from the call stack would be computationally inefficient
%
% Syntax:
%    [res,vars] = pre_norm(classname,S)
%    [res,vars] = pre_norm(classname,S,type)
%    [res,vars] = pre_norm(classname,S,type,mode)
%
% Inputs:
%    classname - class of calling function
%    S - contSet object
%    type - (optional) norm type (default: 2)
%    mode - (optional) mode (default: 'ub')
%
% Outputs:
%    res - true/false whether result has been computed here
%    vars - cell-array whose content depends on value of res:
%           res = true: value of norm
%           res = false: checked input argument arguments incl. default
%                        values for norm function evaluation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-August-2022
% Last update:  23-November-2022 (MW, add classname as input argument)
% Last revision:---

%------------- BEGIN CODE --------------

% check number of input arguments
if nargin < 1
    throw(CORAerror('CORA:notEnoughInputArgs',1));
elseif nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs',4));
end

% set default values: type (2) and mode ('ub')
[type,mode] = setDefaultValues({2,'ub'},varargin);

% check input arguments: two tries depending on value for type
try
    inputArgsCheck({{S,'att',classname},...
                    {type,'att','numeric','scalar'},...
                    {mode,'str',{'exact','ub','ub_convex'}}});
    % have to check value of type
    if ~any(type == [1,2,Inf])
        throw(CORAerror('CORA:wrongValue','second','1, 2, or Inf'));
    end
catch
    % type not numeric (if issue with S/mode, then same issue here; time
    % consumption does not matter since code cannot proceed anyway)
    inputArgsCheck({{S,'att',classname}, ...
                    {type,'str','fro'},...
                    {mode,'str',{'exact','ub','ub_convex'}}});
end

% initialize that no result has been computed
res = false;

% special case: contSet is empty set
if isemptyobject(S)
    res = true;
    vars{1} = -Inf; return
end

% result has to be computed in calling function
vars{1} = S;
vars{2} = type;
vars{3} = mode;

%------------- END OF CODE --------------
