function res = test_partition_intersection()
% test_partition_intersection - unit test for the overapproximating way of 
% finding segments in a partition intersected by a continuous set, i.e. the 
% function intersectingCells
%
% Syntax:
%    res = test_partition_intersection()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Aaron Pereira, Matthias Althoff
% Written:       02-August-2017
% Last update:   02-August-2018 (MA)
%                18-February-2025 (MA)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate three-dimensional field
threeDimField_div = partition({[0 2 3 4 8 10],...
    [-3 -1.5 -1 -0.9 0 0.1 0.2 0.3 1 2 3],...
    [0,0.3,0.6,1]});

% create the same set as interval, polytope, and zonotope inside the
% partitioned field
I = interval([1;0.2;0.1],[4;2.2;0.6]);
P = polytope(I);
Z = zonotope(I);

% return cells intersected by the polytope, interval, and zonotope
intCellsP = intersectingCells(threeDimField_div,P);
intCellsI = intersectingCells(threeDimField_div,I);
intCellsZ = intersectingCells(threeDimField_div,Z);

% check if the same cells are intersected
assert((~any(unique(intCellsP) - unique(intCellsZ)))...
     &&(~any(unique(intCellsI) - unique(intCellsZ))));

% create the same set as interval, polytope, and zonotope partly outside the
% partitioned field
I = interval([1;0.2;-0.1],[4;2.2;0.6]);
P = polytope(I);
Z = zonotope(I);

% return cells intersected by the polytope, interval, and zonotope
intCellsP = intersectingCells(threeDimField_div,P);
intCellsI = intersectingCells(threeDimField_div,I);
intCellsZ = intersectingCells(threeDimField_div,Z);

% check if the same cells are intersected
assert((~any(unique(intCellsP) - unique(intCellsZ)))&&...
       (~any(unique(intCellsI) - unique(intCellsZ))));

% first compute cell segments and then cell indices to check whether the
% conversion from segments to indices works correctly
% polytope
intSSP = intersectingCells(threeDimField_div,P,'subscripts');
intCellsP1 = cellIndices(threeDimField_div,intSSP);
% interval
intSSI = intersectingCells(threeDimField_div,I,'subscripts');
intCellsI1 = cellIndices(threeDimField_div,intSSI);
% zonotope
intSSZ = intersectingCells(threeDimField_div,Z,'subscripts');
intCellsZ1 = cellIndices(threeDimField_div,intSSZ);

% check if the same cells are intersected
assert((~any(unique(intCellsP1)-unique(intCellsP)))&&...
       (~any(unique(intCellsI1)-unique(intCellsI)))&&...
       (~any(unique(intCellsZ1)-unique(intCellsZ))));

% the test is passed if all assertions are correct
res = true;

% ------------------------------ END OF CODE ------------------------------
