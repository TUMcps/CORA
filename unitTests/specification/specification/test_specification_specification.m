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

% assume true
res = true;

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
if length(spec) ~= 1 || ~isequal(spec.set,set) || ~strcmp(spec.type,'unsafeSet')
    res = false;
end

% list of sets
spec = specification(list);
if length(spec) ~= 3 || ~isequal(spec(1).set,list{1}) ...
        || ~isequal(spec(2).set,list{2}) || ~isequal(spec(3).set,list{3}) ...
        || ~strcmp(spec(1).type,'unsafeSet') || ~strcmp(spec(2).type,'unsafeSet') ...
        || ~strcmp(spec(3).type,'unsafeSet')
    res = false;
end

% list of sets with different type
spec = specification(list,'safeSet');
if length(spec) ~= 3 || ~isequal(spec(1).set,list{1}) ...
        || ~isequal(spec(2).set,list{2}) || ~isequal(spec(3).set,list{3}) ...
        || ~strcmp(spec(1).type,'safeSet') || ~strcmp(spec(2).type,'safeSet') ...
        || ~strcmp(spec(3).type,'safeSet')
    res = false;
end

% list of sets with location
spec = specification(list,'safeSet',location_pHA);
if length(spec) ~= 3 || ~isequal(spec(1).set,list{1}) ...
        || ~isequal(spec(2).set,list{2}) || ~isequal(spec(3).set,list{3}) ...
        || ~strcmp(spec(1).type,'safeSet') || ~strcmp(spec(2).type,'safeSet') ...
        || ~strcmp(spec(3).type,'safeSet') || ~isequal(spec(1).location,location_pHA) ...
        || ~isequal(spec(2).location,location_pHA) || ~isequal(spec(3).location,location_pHA)
    res = false;
end

% set with type and time 
spec = specification(set,'invariant',time);
if length(spec) ~= 1 || ~isequal(spec.set,set) || ~strcmp(spec.type,'invariant') ...
        || ~isequal(spec.time,time)
    res = false;
end

% list of sets with type and time 
spec = specification(list,'safeSet',time);
if length(spec) ~= 3 || ~isequal(spec(1).set,list{1}) ...
        || ~isequal(spec(2).set,list{2}) || ~isequal(spec(3).set,list{3}) ...
        || ~strcmp(spec(1).type,'safeSet') || ~strcmp(spec(2).type,'safeSet') ...
        || ~strcmp(spec(3).type,'safeSet') || ~isequal(spec(1).time,time) ...
        || ~isequal(spec(2).time,time) || ~isequal(spec(3).time,time)
    res = false;
end

% list of sets with type, location, and time
spec = specification(list,'safeSet',time,location_pHA);
if length(spec) ~= 3 || ~isequal(spec(1).set,list{1}) ... 
        || ~isequal(spec(2).set,list{2}) || ~isequal(spec(3).set,list{3}) ...
        || ~strcmp(spec(1).type,'safeSet') || ~strcmp(spec(2).type,'safeSet') ...
        || ~strcmp(spec(3).type,'safeSet') || ~isequal(spec(1).time,time) ...
        || ~isequal(spec(3).time,time) || ~isequal(spec(3).time,time) ...
        || ~isequal(spec(1).location,location_pHA) ...
        || ~isequal(spec(2).location,location_pHA) ...
        || ~isequal(spec(3).location,location_pHA)
    res = false;
end

% stl formula
spec = specification(eq);
if ~strcmp(spec.type,'logic')
    res = false;
end
spec = specification(eq,'logic');
if ~strcmp(spec.type,'logic')
    res = false;
end

% function handle
spec = specification(funHan);
if ~isequal(spec.set,funHan) || ~strcmp(spec.type,'custom')
    res = false;
end
spec = specification(funHan,'custom');
if ~isequal(spec.set,funHan) || ~strcmp(spec.type,'custom')
    res = false;
end


% wrong instantiations

% stl formula with wrong type
try
    spec = specification(eq,'safeSet');
    res = false;
end

% stl formula with too many input arguments
try
    spec = specification(eq,'logic',time);
    res = false;
end
% stl formula with too many input arguments
try
    spec = specification(eq,'logic',location_HA);
    res = false;
end

% function handle with wrong type
try
    spec = specification(funHan,'safeSet');
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
