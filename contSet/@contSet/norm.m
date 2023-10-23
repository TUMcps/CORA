function varargout = norm(S,varargin)
% norm - compute the norm of a set
%
% Syntax:
%    [res,x] = norm(S)
%    [res,x] = norm(S,type)
%    [res,x] = norm(S,type,mode)
%
% Inputs:
%    classname - class of calling function
%    S - contSet object
%    type - (optional) norm type (default: 2)
%    mode - (optional) mode (default: 'ub')
%
% Outputs:
%    res - norm
%    x - point where norm is attained
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       17-August-2022
% Last update:   23-November-2022 (MW, add classname as input argument)
% Last revision: 27-March-2023 (MW, restructure relation to subclass)

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin < 1
    throw(CORAerror('CORA:notEnoughInputArgs',1));
elseif nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% set default values: type (2) and mode ('ub')
[type,mode] = setDefaultValues({2,'ub'},varargin);

% check input arguments: two tries depending on value for type
try
    inputArgsCheck({{S,'att','contSet'},...
                    {type,'att','numeric','scalar'},...
                    {mode,'str',{'exact','ub','ub_convex'}}});
    % have to check value of type
    if ~any(type == [1,2,Inf])
        throw(CORAerror('CORA:wrongValue','second','1, 2, or Inf'));
    end
catch
    % type not numeric (if issue with S/mode, then same issue here; time
    % consumption does not matter since code cannot proceed anyway)
    inputArgsCheck({{S,'att','contSet'}, ...
                    {type,'str','fro'},...
                    {mode,'str',{'exact','ub','ub_convex'}}});
end

% call subclass method
varargout = cell(1,max(1,nargout));
try
    [varargout{:}] = norm_(S,type,mode);
catch ME
    % empty set case
    if representsa_(S,'emptySet',eps)
        varargout{1} = -Inf;
        varargout{2} = [];
    else
        rethrow(ME);
    end

end

% ------------------------------ END OF CODE ------------------------------
