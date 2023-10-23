function res = isempty(S)
% isempty - (DEPRECATED -> representsa(S,'emptySet')
%
% Syntax:
%    res = isempty(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: representsa

% Authors:       Mark Wetzlinger
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check call-stack for workspace calls...
st = dbstack("-completenames");
for i=1:length(st)
    if strcmp(st(i).name,'workspacefunc') || strcmp(st(i).name,'datatipinfo')
        res = builtin('isempty',S); % call builtin isempty
        return
    end
end

warning(sprintf(['The function ''isempty'' is deprecated (since CORA 2024) and has been replaced by ''representsa''.\n' ...
    '         The main reason is that the function ''isempty'' is also called implicitly by MATLAB in various circumstances.\n' ...
    '         For contSet classes with the set property ''emptySet'', this would return different results depending on running or debugging.\n' ...
    '         When updating the code, please rename every function call ''isempty(S)'' -> ''representsa(S,''emptySet'')''.\n' ...
    '         The function ''representsa'' also supports comparison to other contSet classes or special representations.\n' ...
    '         Note that the function ''isempty'' for contSet object will be removed in a future release.']));
res = representsa_(S,'emptySet',eps);

% ------------------------------ END OF CODE ------------------------------
