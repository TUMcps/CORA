function res = isequalFunctionHandle(f,g)
% isequalFunctionHandle - checks if two function handles are equal
%
% Syntax:
%    res = isequalFunctionHandle(f,g)
%
% Inputs:
%    f - function handle
%    g - function handle
%
% Outputs:
%    res - true/false
%
% Example:
%    f = @(x,u) [x(1) - u(1); x(1)*x(2)];
%    g = @(x,u) [x(1) - u(1); x(2)*x(1)];
%    isequalFunctionHandle(f,g);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       07-October-2024 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that we have two function handles
narginchk(2,2);
inputArgsCheck({{f,'att','function_handle'}; ...
                {g,'att','function_handle'}});

% get the number of input/output arguments for each function handle
[f_in,f_out] = inputArgsLength(f);
[g_in,g_out] = inputArgsLength(g);

% number of input arguments and their size must match
if any(size(f_in) ~= size(g_in)) || any(f_in ~= g_in)
    res = false;
    return
end

% number of output arguments must match
if f_out ~= g_out
    res = false;
    return
end

% insert symbolic variables into the functions
numArgsIn = numel(f_in);
argsIn = cell(numArgsIn,1);
for i=1:numArgsIn
    varName = ['argsIn_' num2str(i) '_'];
    argsIn{i} = sym(varName,[f_in(i),1],'real');
end
f_sym = f(argsIn{:});
g_sym = g(argsIn{:});

% check symbolic expressions for equality
res = isequal(f_sym,g_sym);

% ------------------------------ END OF CODE ------------------------------
