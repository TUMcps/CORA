function res = test_partition_intersection()
% test_partition_intersection - unit test for the overapproximating way of 
% finding segments in a partition intersected by a continuous set, i.e. the 
% function intersectingCells
%
% Syntax:  
%    res = test_partition_intersection()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Aaron Pereira, Matthias Althoff
% Written:      02-August-2017
% Last update:  02-August-2018 (MA)
% Last revision:---

%------------- BEGIN CODE --------------

threeDimField_div=partition({[0 2 3 4 8 10],[-3 -1.5 -1 -0.9 0 0.1 0.2 0.3 1 2 3],[0,0.3,0.6,1]});

P = mptPolytope([2 0.2 0.3;4 2.2 0.6;1 1.2 0.5; 1 1.2 0.1]);
I = interval([1;0.2;0.1],[4;2.2;0.6]);
Z = zonotope([2.5 1.2 0.35;1.5 0 0;0 1 0;0 0 0.25]');

intCellsP = intersectingCells(threeDimField_div,P);
intCellsI = intersectingCells(threeDimField_div,I);
intCellsZ = intersectingCells(threeDimField_div,Z);

if (length(intCellsP) == length(intCellsI))&&(length(intCellsZ) == length(intCellsI))&&(length(intCellsP) == length(intCellsZ))
    res1 = (~any(unique(intCellsP) - unique(intCellsZ)))&&(~any(unique(intCellsI) - unique(intCellsZ)));
else
    res1 = 0;
end

% when slightly outside!

P = mptPolytope([2 0.2 0.3;4 2.2 0.6;1 1.2 0.5; 1 1.2 -0.1]);
I = interval([1;0.2;-0.1],[4;2.2;0.6]);
Z = zonotope([2.5 1.2 0.25;1.5 0 0;0 1 0;0 0 0.35]');

intCellsP = intersectingCells(threeDimField_div,P);
intCellsI = intersectingCells(threeDimField_div,I);
intCellsZ = intersectingCells(threeDimField_div,Z);

if (length(intCellsP) == length(intCellsI))&&(length(intCellsZ) == length(intCellsI))&&(length(intCellsP) == length(intCellsZ))
    res2 = (~any(unique(intCellsP) - unique(intCellsZ)))&&(~any(unique(intCellsI) - unique(intCellsZ)));
else
    res2 = 0;
end

intSSP = intersectingCells(threeDimField_div,P,'subscripts');
intCellsP1 = cellIndices(threeDimField_div,intSSP);
intSSI = intersectingCells(threeDimField_div,I,'subscripts');
intCellsI1 = cellIndices(threeDimField_div,intSSI);
intSSZ = intersectingCells(threeDimField_div,Z,'subscripts');
intCellsZ1 = cellIndices(threeDimField_div,intSSZ);

if (length(intCellsP1)==length(intCellsP))&&(length(intCellsI1)==length(intCellsI))&&(length(intCellsZ1)==length(intCellsZ))
    res3 = (~any(unique(intCellsP1)-unique(intCellsP)))&&(~any(unique(intCellsI1)-unique(intCellsI)))&&(~any(unique(intCellsZ1)-unique(intCellsZ)));
else
    res3 = 0;
end

res = res1&&res2&&res3;

%------------- END OF CODE --------------