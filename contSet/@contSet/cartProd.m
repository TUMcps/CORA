function res = cartProd(S1,S2,varargin)
% cartProd - computes the Cartesian product of two sets
%
% Description:
%    computes the set { [s_1 s_2 ]^T | s_1 \in \mathcal{S}_1, s_2 \in \mathcal{S}_2 }.
%
% Syntax:
%    res = cartProd(S1,S2)
%    res = cartProd(S1,S2,type)
%
% Inputs:
%    S1,S2 - contSet object
%    type - type of computation ('exact','inner','outer')
%
% Outputs:
%    res - Cartesian product
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       18-August-2022
% Last update:   23-November-2022 (MW, add classname as input argument)
%                03-January-2023 (MW, fix bug regarding reordering of args)
% Last revision: 27-March-2023 (MW, restructure relation to subclass)

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs',2));
elseif nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
if isa(S1,'ellipsoid') || isa(S2,'ellipsoid')
    type = setDefaultValues({'outer'},varargin);
else
    type = setDefaultValues({'exact'},varargin);
end

% check input arguments: two versions as order of input arguments matters
try
    inputArgsCheck({{S1,'att',{'contSet','numeric'}};
                    {S2,'att',{'contSet','numeric'},'vector'};
                    {type,'str',{'outer','inner','exact'}}});
catch
    inputArgsCheck({{S1,'att',{'contSet','numeric'},'vector'};
                    {S2,'att',{'contSet','numeric'}};
                    {type,'str',{'outer','inner','exact'}}});
end

% call subclass method
try
    res = cartProd_(S1,S2,type);
catch ME
    % Cartesian products with empty sets are currently not supported,
    % because we cannot concatenate empty vectors with filled vectors
    if representsa_(S1,'emptySet',1e-8) || representsa_(S2,'emptySet',1e-8)
        throw(CORAerror('CORA:notSupported',...
            "Cartesian products with empty sets are not supported."));
    else
        rethrow(ME);
    end

end

% ------------------------------ END OF CODE ------------------------------
