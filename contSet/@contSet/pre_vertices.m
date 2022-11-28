function [res,vars] = pre_vertices(classname,S,varargin)
% pre_vertices - pre-processing for contSet/vertices functions including
%    1. setting of default values
%    2. checking of input arguments
%    3. computation of special cases (e.g., empty set)
%    the class name of the calling function is handed as an input argument
%    as a read-out from the call stack would be computationally inefficient
%
% Syntax:
%    [res,vars] = pre_vertices(classname,S)
%    [res,vars] = pre_vertices(classname,S,method)
%
% Inputs:
%    classname - class of calling function
%    S - contSet object
%    method - method for computation of vertices
%
% Outputs:
%    res - true/false whether result has been computed here
%    vars - cell-array whose content depends on value of res:
%           res = true: value of vertices
%           res = false: checked input argument arguments incl. default
%                        values for vertices function evaluation
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
elseif nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
method = setDefaultValues({'convHull'},varargin{:}); 

% check input arguments
inputArgsCheck({{S,'att',classname};
                {method,'str',{'convHull','iterate','polytope'}}});

% empty set case
if isemptyobject(S)
    res = true;
    vars{1} = [];
    return
end

% result has to be computed in calling function
res = false;
vars = {S,method};

%------------- END OF CODE --------------
