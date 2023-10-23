function res = test_specification_add
% test_specification_add - unit test for add
%
% Syntax:
%    res = test_specification_add
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
set3 = capsule([1;1],[1;1],0.5);
list = {set,set2,set3};

% init time
time = interval(2,4);

% init location
location_HA = [1,2];
location_pHA = {[1,2], 2, [2,3]};

% init stl formula
% x = stl('x',2);
% eq = until(x(1) <= 5,x(2) > 3 & x(1) <= 2, interval(0.1,0.2));

% init function handle
% funHan = @(x) x(1)^2;

% one set
spec{1} = specification(set);
% list of sets
spec{2} = specification(list);
% list of sets with different type
spec{3} = specification(list,'safeSet');
% list of sets with location
spec{4} = specification(list,'safeSet',location_HA);
% set with type and time
spec{5} = specification(set,'invariant',time);
% list of sets with type and time 
spec{6} = specification(list,'safeSet',time);
% list of sets with type, location, and time
spec{7} = specification(list,'safeSet',time,location_pHA);
% stl formula
%spec{8} = specification(eq);
% function handle
%spec{9} = specification(funHan);


% add all combinations of specifications
comb = combinator(length(spec),2,'c');
for i=1:size(comb,1)
    % read out specification objects
    spec1 = spec{comb(i,1)};
    spec2 = spec{comb(i,2)};

    % join specification objects
    spec_add = add(spec1,spec2);

    % compare specification objects
    for j=1:length(spec1)+length(spec2)
        if j <= length(spec1)
            % from first part
            if ~isequal(spec_add(j),spec1(j))
                res = false; return
            end
        
        else
            % from second part
            if ~isequal(spec_add(j),spec2(j-length(spec1)))
                res = false; return
            end

        end
    end
end

% ------------------------------ END OF CODE ------------------------------
