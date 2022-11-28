function [res,vars] = pre_cartProd(classname,S1,S2,varargin)
% pre_cartProd - pre-processing for contSet/cartProd functions including
%    1. setting of default values
%    2. checking of input arguments
%    3. computation of special cases (e.g., empty set)
%    the class name of the calling function is handed as an input argument
%    as a read-out from the call stack would be computationally inefficient
%
% Syntax:
%    [res,vars] = pre_cartProd(classname,S1,S2)
%    [res,vars] = pre_cartProd(classname,S1,S2,type)
%
% Inputs:
%    classname - class of calling function
%    S1,S2 - contSet object
%    type - type of computation ('exact','inner','outer')
%
% Outputs:
%    res - true/false whether result has been computed here
%    vars - cell-array whose content depends on value of res:
%           res = true: value of cartProd
%           res = false: checked input argument arguments incl. default
%                        values for cartProd function evaluation
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
if nargin < 3
    throw(CORAerror('CORA:notEnoughInputArgs',3));
elseif nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs',4));
end

% right order of objects
[S1,S2] = findClassArg(S1,S2,classname);

% parse input arguments
if strcmp(classname,'ellipsoid')
    type = setDefaultValues({'outer'},varargin{:});
else
    type = setDefaultValues({'exact'},varargin{:});
end

% check input arguments
inputArgsCheck({{S1,'att',classname};
                {S2,'att',{'contSet','numeric'},'vector'};
                {type,'str',{'outer','inner','exact'}}});


% empty set cases (current handling not entirely mathematically correct)
if isemptyobject(S1)
    res = true;
    vars{1} = S2;
    return
elseif (isnumeric(S2) && isempty(S2)) || (isa(S2,'contSet') && isemptyobject(S2))
    res = true;
    vars{1} = S1;
    return
end

% % Cartesian products with empty sets are not supported
% % if isemptyobject(S1) || ...
% %     (isnumeric(S2) && isempty(S2)) || (isa(S2,'contSet') && isemptyobject(S2))
% %     throw(CORAerror('CORA:notSupported',...
% %         "Cartesian products with empty sets are not supported."));
% % end

% result has to be computed in calling function
res = false;
vars = {S1,S2,type};

%------------- END OF CODE --------------
