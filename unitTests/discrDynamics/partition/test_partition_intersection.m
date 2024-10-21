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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

threeDimField_div=partition({[0 2 3 4 8 10],[-3 -1.5 -1 -0.9 0 0.1 0.2 0.3 1 2 3],[0,0.3,0.6,1]});

P = polytope([2 0.2 0.3;4 2.2 0.6;1 1.2 0.5; 1 1.2 0.1]');
I = interval([1;0.2;0.1],[4;2.2;0.6]);
Z = zonotope([2.5 1.2 0.35;1.5 0 0;0 1 0;0 0 0.25]');

intCellsP = intersectingCells(threeDimField_div,P);
intCellsI = intersectingCells(threeDimField_div,I);
intCellsZ = intersectingCells(threeDimField_div,Z);

assert((length(intCellsP) == length(intCellsI))&&(length(intCellsZ) == length(intCellsI))&&(length(intCellsP) == length(intCellsZ)))
assert((~any(unique(intCellsP) - unique(intCellsZ)))&&(~any(unique(intCellsI) - unique(intCellsZ))));

% when slightly outside!

P = polytope([2 0.2 0.3;4 2.2 0.6;1 1.2 0.5; 1 1.2 -0.1]');
I = interval([1;0.2;-0.1],[4;2.2;0.6]);
Z = zonotope([2.5 1.2 0.25;1.5 0 0;0 1 0;0 0 0.35]');

intCellsP = intersectingCells(threeDimField_div,P);
intCellsI = intersectingCells(threeDimField_div,I);
intCellsZ = intersectingCells(threeDimField_div,Z);

assert((length(intCellsP) == length(intCellsI))&&(length(intCellsZ) == length(intCellsI))&&(length(intCellsP) == length(intCellsZ)))
assert((~any(unique(intCellsP) - unique(intCellsZ)))&&(~any(unique(intCellsI) - unique(intCellsZ))));

intSSP = intersectingCells(threeDimField_div,P,'subscripts');
intCellsP1 = cellIndices(threeDimField_div,intSSP);
intSSI = intersectingCells(threeDimField_div,I,'subscripts');
intCellsI1 = cellIndices(threeDimField_div,intSSI);
intSSZ = intersectingCells(threeDimField_div,Z,'subscripts');
intCellsZ1 = cellIndices(threeDimField_div,intSSZ);

assert((length(intCellsP1)==length(intCellsP))&&(length(intCellsI1)==length(intCellsI))&&(length(intCellsZ1)==length(intCellsZ)))
assert((~any(unique(intCellsP1)-unique(intCellsP)))&&(~any(unique(intCellsI1)-unique(intCellsI)))&&(~any(unique(intCellsZ1)-unique(intCellsZ))));

res = true;

% ------------------------------ END OF CODE ------------------------------
