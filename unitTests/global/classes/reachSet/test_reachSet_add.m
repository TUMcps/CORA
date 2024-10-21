function res = test_reachSet_add
% test_reachSet_add - unit test function for add
%
% Syntax:
%    res = test_reachSet_add()
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
% Written:       06-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% instantiate nx1 reachSet objects with different parents
timePoint.set{1} = zonotope(zeros(2,1),eye(2));
timePoint.time{1} = 0;

parent1 = [0; 1; 1];
parent2 = [0; 1; 0; 2; 2];
R1 = []; R2 = [];
lengthR1 = 3;
for i=1:lengthR1
    R1 = [R1; reachSet(timePoint,[],parent1(i))];
end
for i=1:5
    R2 = [R2; reachSet(timePoint,[],parent2(i))];
end


% add empty reachSet object
assert(isequal(R1,add(R1,reachSet())))
assert(isequal(R1,add(reachSet(),R1)))


% add reachSet objects
R = add(R1,R2);

% parent property of R1 remains the same
for i=1:length(R1)
    assertLoop(R(i).parent == parent1(i),i)
end
% parent property of R2 incremented by length of R1 only if non-zero
for i=1:length(R2)
    if R2(i).parent == 0
        assertLoop(R(i+lengthR1).parent == 0,i);
    else
        assertLoop(R(i+lengthR1).parent == R2(i).parent + lengthR1,i)
    end
end


% increment by non-zero parent
shift = 3;
R = add(R1,R2,shift);

% parent property of R1 remains the same
for i=1:length(R1)
    assertLoop(R(i).parent == parent1(i),i)
end
% parent property of R2 incremented by length of R1 only if non-zero
for i=1:length(R2)
    if R2(i).parent == 0
        assertLoop(R(i+lengthR1).parent == shift,i)
    else
        assertLoop(R(i+lengthR1).parent == R2(i).parent + lengthR1,i)
    end
end

% ------------------------------ END OF CODE ------------------------------
