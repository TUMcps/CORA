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

CORAwarning('CORA:deprecated','function','contSet/isempty','CORA v2024', ...
    'When updating the code, please replace every function call ''isempty(S)'' with ''representsa(S,''emptySet'')''.', ...
    ['The main reason is that the function ''isempty'' is also called implicitly by MATLAB in various circumstances.\n' ...
    'For contSet classes with the set property ''emptySet'', this would return different results depending on running or debugging.'])
res = representsa_(S,'emptySet',eps);

% ------------------------------ END OF CODE ------------------------------
