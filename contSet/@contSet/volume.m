function res = volume(S,varargin)
% volume - computes the volume of a set
%
% Syntax:
%    res = volume(S)
%    res = volume(S,method)
%    res = volume(S,method,order)
%
% Inputs:
%    S - contSet object
%    method - (optional) method for evaluation
%    order - (optional) zonotope order
%
% Outputs:
%    res - volume
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
if nargin < 1
    throw(CORAerror('CORA:notEnoughInputArgs',1));
elseif nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% check input arguments
inputArgsCheck({{S,'att','contSet'}});

% only zonotopes: set default values method ('exact') and order (5)
[method,order] = setDefaultValues({'exact',5},varargin);
% check input arguments
inputArgsCheck({{method,'str',{'exact','reduce','alamo'}},...
                {order,'att','numeric',{'integer','positive'}}});

% call subclass method
try
    res = volume_(S,method,order);
catch ME
    % empty set case
    if representsa_(S,'emptySet',eps)
        res = 0;
    else
        rethrow(ME);
    end
end

% ------------------------------ END OF CODE ------------------------------
