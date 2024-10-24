function res = test_specification_project
% test_specification_project - unit test for project
%
% Syntax:
%    res = test_specification_project
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
set = zonotope([0;0;0],[1,-0.7,0.3;0.2,1,-1;0,0.2,-0.2]);
set2 = interval([-1;-2;0],[2;3;1]);
list = {set,set2};

% project sets
projDim = [1,3];
set_ = project(set,projDim);
set2_ = project(set2,projDim);
list_ = {set_,set2_};

% init time
time = interval(2,4);

% init location
location_HA = [1,2];

% init stl formula
x = stl('x',2);
eq = until(x(1) <= 5,x(2) > 3 & x(1) <= 2, interval(0.1,0.2));

% init function handle
funHan = @(x) x(1)^2 - x(2)*x(3);

% single set
spec = specification(set,'unsafeSet');
spec_ = project(spec,projDim);
% compare set
assert(isequal(spec_.set,set_));

% one unsafe set
spec = specification(set,'safeSet');
spec_ = project(spec,projDim);
% compare set
assert(isequal(spec_.set,set_));

% list of unsafe sets
spec = specification(list);
spec_ = project(spec,projDim);
% check all sets
for i=1:length(spec_)
    assert(isequal(spec_(i).set,list_{i}));
end

% list of invariants
spec = specification(list,'invariant');
spec_ = project(spec,projDim);
% check all sets
for i=1:length(spec_)
    % type has to stay the same, check projection
    assert(strcmp(spec_(i).type,'invariant'))
    assert(isequal(spec_(i).set,list_{i}));
end

% list of sets with time and location
spec = specification(list,'safeSet',time,location_HA);
spec_ = project(spec,projDim);
for i=1:length(spec_)
    % type, time, and location has to stay the same, check projection
    assert(strcmp(spec_(i).type,'safeSet'))
    assert(isequal(spec_(i).time,spec(i).time))
    assert(all(spec_(i).location == spec(i).location))
    assert(isequal(spec_(i).set,list_{i}));
end

% stl formula
spec = specification(eq,'logic');
% 'logic' not supported
assertThrowsAs(@project,'CORA:notSupported',spec,projDim);

% function handle
spec = specification(funHan,'custom');
% 'custom' not supported
assertThrowsAs(@project,'CORA:notSupported',spec,projDim);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
