function res = verify(obj,params,options,spec)
% verify - verifies a hybrid system against the given specifications
%
% Syntax:
%    R = verify(obj,params,options,spec)
%
% Inputs:
%    obj - hybridDynamics object
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%    spec - object of class specification or stl
%
% Outputs:
%    res - true/false whether specifications are satisfied
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics/verify

% Authors:       Tobias Ladner
% Written:       19-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. check number of inputs
if nargin < 4
    throw(CORAerror('CORA:notEnoughInputArgs',4))
elseif nargin > 4
    throw(CORAerror("CORA:tooManyInputArgs", 4))
end

% 2. validate inputs
inputArgsCheck({ ...
    {obj, 'att', 'hybridDynamics'}; ...
    {params, 'att', 'struct'}; ...
    {options, 'att', 'struct'}; ...
    {spec, 'att', {'specification','stl'}}; ...
});

% 3. call verify function depending on given specification
switch class(spec)
    case 'stl'
        res = verifySTL_seidl(obj,params,options,spec);

    otherwise
        throw(CORAerror('CORA:noops',obj,params,options,spec));
end


end

% ------------------------------ END OF CODE ------------------------------
