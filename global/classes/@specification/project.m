function spec = project(spec,dims)
% project - projects the set of a specification onto a subspace
%
% Syntax:  
%    spec = project(spec,dims)
%
% Inputs:
%    spec - specification object
%    dims - dimensions for projection
%
% Outputs:
%    spec - projected specification object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Mark Wetzlinger
% Written:      26-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% project set
switch spec.type
    case {'invariant','safeSet','unsafeSet'}
        % loop over array of specifications
        for i=1:length(spec)
            spec(i).set = project(spec(i).set,dims);
        end
    case 'custom'
        throw(CORAerror('CORA:notSupported',...
            "Projection a specifications of type 'custom' is not yet supported."));
end

%------------- END OF CODE --------------