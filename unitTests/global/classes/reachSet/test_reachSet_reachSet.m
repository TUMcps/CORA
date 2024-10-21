function res = test_reachSet_reachSet
% test_reachSet_reachSet - unit test function for constructor
%
% Syntax:
%    res = test_reachSet_reachSet()
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
% Written:       11-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty reachSet
R = reachSet();

% initialize some sets for instantiations
Z = zonotope([1; 1],[1 0 3; -2 1 4]);

timePoint.set{1,1} = Z;
timePoint.time{1,1} = 0;
for i=1:5
    timePoint.set{i+1,1} = Z;
    timePoint.time{i+1,1} = i;
    timeInt.set{i,1} = Z;
    timeInt.time{i,1} = interval((i-1),i);
end
parent = 1;
loc = 1;

% non-matching lengths
timePoint_.set{1,1} = Z;
timePoint_.time{1,1} = 0;
for i=1:5
    timePoint_.set{i+1,1} = Z;
    timeInt_.time{i,1} = interval((i-1),i);
    if i <= 4
        timePoint_.time{i+1,1} = i;
        timeInt_.set{i,1} = Z;
    end
end
parent_ = -1;
loc_ = -1;

% correct instantiations according to constructor (all with correct length)
R = reachSet(timePoint);
R = reachSet(timePoint,parent);
R = reachSet(timePoint,parent,loc);
R = reachSet(timePoint,timeInt);
R = reachSet(timePoint,timeInt,parent);
R = reachSet(timePoint,timeInt,parent,loc);
% non-computed sets
R = reachSet([],timeInt);
R = reachSet([],timeInt,parent);
R = reachSet(timePoint,[]);
R = reachSet(timePoint,[],parent);

% wrong instantiations: non-matching lengths, negative values, too many
% input arguments
assertThrowsAs(@reachSet,'CORA:wrongInputInConstructor',timePoint_);
assertThrowsAs(@reachSet,'CORA:wrongInputInConstructor',timePoint,timeInt_);
assertThrowsAs(@reachSet,'MATLAB:validators:mustBeNonnegative',timePoint,timeInt,parent_);
assertThrowsAs(@reachSet,'MATLAB:validators:mustBeNonnegative',timePoint,timeInt,parent,loc_);
assertThrowsAs(@reachSet,'CORA:numInputArgsConstructor',timePoint,timeInt,parent,loc,loc);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
