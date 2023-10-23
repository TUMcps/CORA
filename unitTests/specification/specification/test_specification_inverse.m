function res = test_specification_inverse
% test_specification_inverse - unit test for inverse
%
% Syntax:
%    res = test_specification_inverse
%
% Inputs:
%    -
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
% Written:       30-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% init sets
set = zonotope([0;0],[1,-0.7;0.2,1]);
set2 = interval([-1;-2],[2;3]);
list = {set,set2};

% init time
time = interval(2,4);

% init location
location_HA = [1,2];

% init stl formula
x = stl('x',2);
eq = until(x(1) <= 5,x(2) > 3 & x(1) <= 2, interval(0.1,0.2));

% init function handle
funHan = @(x) x(1)^2;

% one unsafe set
spec = specification(set,'unsafeSet');
spec_ = inverse(spec);
% single 'unsafeSet' becomes 'safeSet'
res(end+1,1) = strcmp(spec_.type,'safeSet');

% one unsafe set
spec = specification(set,'safeSet');
spec_ = inverse(spec);
% 'safeSet' becomes 'unsafeSet'
res(end+1,1) = strcmp(spec_.type,'unsafeSet');

% list of unsafe sets
spec = specification(list);
spec_ = inverse(spec);
% list of unsafe sets remain unsafe
res(end+1,1) = all(arrayfun(@(x) strcmp(x.type,'unsafeSet'),spec_,...
    'UniformOutput',true));

% list of safe sets
spec = specification(list,'safeSet');
spec_ = inverse(spec);
% all 'safeSet' become a single 'unsafeSet'
res(end+1,1) = strcmp(spec_.type,'unsafeSet') && length(spec_) == 1;

% list of sets with location
spec = specification(list,'safeSet',location_HA);
spec_ = inverse(spec);
% all 'safeSet' become a single 'unsafeSet'
res(end+1,1) = strcmp(spec_.type,'unsafeSet') && length(spec_) == 1;


% combine results
res = all(res);

% list of sets with type and time 
spec = specification(list,'safeSet',time);
try
    % time not supported
    spec_ = inverse(spec);
    res = false; return
end

% set with type and time
spec = specification(set,'invariant');
try
    % 'invariant' not supported
    spec_ = inverse(spec);
    res = false; return
end

% stl formula
spec = specification(eq,'logic');
try
    % 'logic' not supported
    spec_ = inverse(spec);
    res = false; return
end

% function handle
spec = specification(funHan,'custom');
try
    % 'custom' not supported
    spec_ = inverse(spec);
    res = false; return
end

% ------------------------------ END OF CODE ------------------------------
