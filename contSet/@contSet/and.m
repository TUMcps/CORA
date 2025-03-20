function res = and(S1,S2,varargin)
% and - overloads '&' operator, computes the intersection of two sets
%
% Description:
%    computes the set { s | s \in \mathcal{S}_1, s \in \mathcal{S}_2 }
%
% Syntax:
%    res = S1 & S2
%    res = and(S1,S2)
%    res = and(S1,S2,type)
%
% Inputs:
%    S1,S2 - contSet object
%    type - type of computation ('exact','inner','outer')
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       18-August-2022
% Last update:   23-November-2022 (MW, add classname as input argument)
% Last revision: 27-March-2023 (MW, restructure relation to subclass)

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments (note: andAveraging submethod of
% zonotope/and has up to 6 input arguments, leave that to calling and_)
narginchk(2,3);

% check input arguments
inputArgsCheck({{S1,'att',{'contSet','numeric'}};
                {S2,'att',{'contSet','numeric','cell'}}});

% order input arguments according to their precendence
[S1,S2] = reorder(S1,S2);

if isa(S1,'ellipsoid')
    % parse input arguments
    type = setDefaultValues({'outer'},varargin);

    % check additional input arguments
    inputArgsCheck({{type,'str',{'inner','outer'}}});

elseif isa(S1,'zonotope')
    % parse input arguments
    type = setDefaultValues({'conZonotope'},varargin);

    % check additional input arguments
    inputArgsCheck({{type,'str',{'conZonotope','averaging'}}});
else
    type = setDefaultValues({'exact'},varargin);
end

% check dimension mismatch
equalDimCheck(S1,S2);


% call subclass method
try
    res = and_(S1,S2,type);
catch ME
    if isa(S1,'contSet') && representsa_(S1,'emptySet',eps,'linearize',0,1)
        res = S1;
    elseif isnumeric(S1) && isempty(S1)
        res = [];
    elseif isa(S2,'contSet') && representsa_(S2,'emptySet',eps,'linearize',0,1)
        res = S2;
    elseif isnumeric(S2) && isempty(S2)
        res = [];
    else
        rethrow(ME);
    end

end

% ------------------------------ END OF CODE ------------------------------
