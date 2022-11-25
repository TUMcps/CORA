function res = test_partition_center_functions()
% test_partition_center_functions - unit test for the functions 
% intersectingCells and cellCenter
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

oneDimField=partition([0,10],5);
threeDimField=partition([0,10; -3,3; 0,1],[5;10;3]);
threeDimField_odd=partition([0,10; -3,3; 0,1],[5;11;3]);
threeDimField_div=partition({[0 2 3 4 8 10],[-3 -1.5 -1 -0.9 0 0.1 0.2 0.3 1 2 3],[0,0.3,0.6,1]});

% obtain center cells of partitions
c1 = intersectingCells(oneDimField,sum(oneDimField.intervals,2)/2,'indices');
c2 = intersectingCells(threeDimField,sum(threeDimField.intervals,2)/2,'indices');
c3 = intersectingCells(threeDimField_odd,sum(threeDimField_odd.intervals,2)/2,'indices');
c4 = intersectingCells(threeDimField_div,sum(threeDimField_div.intervals,2)/2,'indices');

% return centers of cells who contain the center of the partition
S1 = cellCenter(oneDimField,c1);
S2 = cellCenter(threeDimField,c2);
S3 = cellCenter(threeDimField_odd,c3);
S4 = cellCenter(threeDimField_div,c4);

res1 = (S1{1}-5)<1e-15;
res2 = norm(S2{1}-[5.0000  -0.3000  0.5000]')<1e-15;
res3 = norm(S2{2}-[5.0000  0.3000  0.5000]')<1e-15;
res4 = norm(S3{1}-[5.0000  0  0.5000]')<1e-15;
res5 = norm(S4{1}-[6.0000  -0.4500   0.4500]')<1e-15;
res6 = norm(S4{2}-[6.0000 0.0500 0.4500]')<1e-15;


res = res1&&res2&&res3&&res4&&res5&&res6;

%------------- END OF CODE --------------