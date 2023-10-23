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
% Example:
%    Z = zonotope([1;-1;0],[1 3 -2; 0 -1 1; 1 2 0]);
%    spec = specification(Z,'safeSet');
%    spec_ = project(spec,[1,3]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       26-June-2022
% Last update:   30-April-2023 (MW, bug fix for arrays, add 'logic')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% loop over array of specifications
for i=1:length(spec)

    % check type
    switch spec(i).type
        case {'invariant','safeSet','unsafeSet'}
            % project set
            spec(i).set = project(spec(i).set,dims);

        case {'custom','logic'}
            throw(CORAerror('CORA:notSupported',...
                "Projection of a specification of types 'custom' or " + ...
                "'logic' is not yet supported."));
    end
end

% ------------------------------ END OF CODE ------------------------------
