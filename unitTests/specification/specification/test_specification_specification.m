function res = test_specification_specification
% test_specification_specification - unit test for constuctor of the class
%    specification
%
% Syntax:
%    res = test_specification_specification
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
% Written:       27-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init sets
set = zonotope([0;0],[1,-0.7;0.2,1]);
set2 = interval([-1;-2],[2;3]);
set3 = capsule([1;1],[1;1],0.5);
list = {set,set2,set3};

% init time
time = interval(2,4);

% init location
location_HA = [1,2];
location_pHA = {[1,2], 2, [2,3]};

% init stl formula
x = stl('x',2);
eq = until(x(1) <= 5,x(2) > 3 & x(1) <= 2, interval(0.1,0.2));

% init function handle
funHan = @(x) x(1)^2;


% instantiate specification objects

% only set
spec = specification(set);
assert(length(spec) == 1)
assert(isequal(spec.set,set))
assert(strcmp(spec.type,'unsafeSet'))

% list of sets
spec = specification(list);
assert(length(spec) == 3)
assert(isequal(spec(1).set,list{1}))
assert(isequal(spec(2).set,list{2}))
assert(isequal(spec(3).set,list{3}))
assert(strcmp(spec(1).type,'unsafeSet'))
assert(strcmp(spec(2).type,'unsafeSet'))
assert(strcmp(spec(3).type,'unsafeSet'))

% list of sets with different type
spec = specification(list,'safeSet');
assert(length(spec) == 3)
assert(isequal(spec(1).set,list{1}))
assert(isequal(spec(2).set,list{2}))
assert(isequal(spec(3).set,list{3}))
assert(strcmp(spec(1).type,'safeSet'))
assert(strcmp(spec(2).type,'safeSet'))
assert(strcmp(spec(3).type,'safeSet'))

% list of sets with location
spec = specification(list,'safeSet',location_pHA);
assert(length(spec) == 3)
assert(isequal(spec(1).set,list{1}))
assert(isequal(spec(2).set,list{2}))
assert(isequal(spec(3).set,list{3}))
assert(strcmp(spec(1).type,'safeSet'))
assert(strcmp(spec(2).type,'safeSet'))
assert(strcmp(spec(3).type,'safeSet'))
assert(isequal(spec(1).location,location_pHA))
assert(isequal(spec(2).location,location_pHA))
assert(isequal(spec(3).location,location_pHA))

% set with type and time 
spec = specification(set,'invariant',time);
assert(length(spec) == 1)
assert(isequal(spec.set,set))
assert(strcmp(spec.type,'invariant'))
assert(isequal(spec.time,time))

% list of sets with type and time 
spec = specification(list,'safeSet',time);
assert(length(spec) == 3)
assert(isequal(spec(1).set,list{1}))
assert(isequal(spec(2).set,list{2}))
assert(isequal(spec(3).set,list{3}))
assert(strcmp(spec(1).type,'safeSet'))
assert(strcmp(spec(2).type,'safeSet'))
assert(strcmp(spec(3).type,'safeSet'))
assert(isequal(spec(1).time,time))
assert(isequal(spec(2).time,time))
assert(isequal(spec(3).time,time))

% list of sets with type, location, and time
spec = specification(list,'safeSet',time,location_pHA);
assert(length(spec) == 3)
assert(isequal(spec(1).set,list{1}))
assert(isequal(spec(2).set,list{2}))
assert(isequal(spec(3).set,list{3}))
assert(strcmp(spec(1).type,'safeSet'))
assert(strcmp(spec(2).type,'safeSet'))
assert(strcmp(spec(3).type,'safeSet'))
assert(isequal(spec(1).time,time))
assert(isequal(spec(3).time,time))
assert(isequal(spec(3).time,time))
assert(isequal(spec(1).location,location_pHA))
assert(isequal(spec(2).location,location_pHA))
assert(isequal(spec(3).location,location_pHA))


% stl formula
spec = specification(eq);
assert(strcmp(spec.type,'logic'));

spec = specification(eq,'logic');
assert(strcmp(spec.type,'logic'));


% function handle
spec = specification(funHan);
assert(isequal(spec.set,funHan));
assert(strcmp(spec.type,'custom'));

spec = specification(funHan,'custom');
assert(isequal(spec.set,funHan));
assert(strcmp(spec.type,'custom'));


% wrong instantiations

% stl formula with wrong type
assertThrowsAs(@specification,'CORA:wrongInputInConstructor',eq,'safeSet');

% stl formula with too many input arguments
assertThrowsAs(@specification,'CORA:notSupported',eq,'logic',time);

% stl formula with too many input arguments
assertThrowsAs(@specification,'CORA:notSupported',eq,'logic',location_HA);

% function handle with wrong type
assertThrowsAs(@specification,'CORA:wrongInputInConstructor',funHan,'safeSet');


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
