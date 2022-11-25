function res = test_partition_exactIntersection()
% test_partition_exactIntersection - unit test for the function 
% exactIntersectingCells
%
% Syntax:  
%    res = test_partition_exactIntersection()
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

twoDimField=partition([0,10; -3,3],[5;10]);
threeDimField=partition([0,10; -3,3; 0,1],[5;10;3]);
threeDimField_div=partition({[0 2 3 4 8 10],[-3 -1.5 -1 -0.9 0 0.1 0.2 0.3 1 2 3],[0,0.3,0.6,1]});

P2 = mptPolytope([2 0;4 2;1 1]);
[iS1,percentages1] = exactIntersectingCells(twoDimField,P2);
% figure;
% plot(twoDimField,iS1(2:end)');
% hold on
% plot(P2);

P2 = mptPolytope([2 0;4 2;-1 1]);
[iS2,percentages2] = exactIntersectingCells(twoDimField,P2);
% figure;
% plot(twoDimField,iS2(2:end)');
% hold on
% plot(P2);

P = mptPolytope([2 0 0.3;4 2 0.6;1 1 0.5; 1 1 0.1]);
[iS3,percentages3] = exactIntersectingCells(threeDimField,P);
% figure;
% plot(threeDimField,iS3(2:end)');
% hold on
% plot(P);

P = mptPolytope([2 0 0.3;4 2 0.6;1 1 0.5; 1 1 -0.1]);
[iS4,percentages4] = exactIntersectingCells(threeDimField,P);
% figure;
% plot(threeDimField,iS4(2:end)');
% hold on
% plot(P);

P = mptPolytope([2 0 0.3;4 2 0.6;1 1 0.5; 1 1 0.1]);
[iS5,percentages5] = exactIntersectingCells(threeDimField_div,P);
% figure;
% plot(threeDimField_div,iS5(2:end)');
% hold on
% plot(P);


P = mptPolytope([2 0 0.3;4 2 0.6;1 1 0.5; 1 1 -0.1]);
[iS6,percentages6] = exactIntersectingCells(threeDimField_div,P);
% figure;
% plot(threeDimField_div,iS6(2:end)');
% hold on
% plot(P);

% iS1_gt = iS1;
% iS2_gt = iS2;
% iS3_gt = iS3;
% iS4_gt = iS4;
% iS5_gt = iS5;
% iS6_gt = iS6;
% percentages1_gt = percentages1;
% percentages2_gt = percentages2;
% percentages3_gt = percentages3;
% percentages4_gt = percentages4;
% percentages5_gt = percentages5;
% percentages6_gt = percentages6;
% 
% save('intersection_ground_truth.mat','iS1_gt','iS2_gt','iS3_gt','iS4_gt','iS5_gt','iS6_gt','percentages1_gt','percentages2_gt','percentages3_gt','percentages4_gt','percentages5_gt','percentages6_gt')

load intersection_ground_truth

accuracy = 1e-14;

res_partial(1) = sum(iS1_gt-iS1) < accuracy ;
res_partial(2) = sum(iS2_gt-iS2) < accuracy ;
res_partial(3) = sum(iS3_gt-iS3) < accuracy ;
res_partial(4) = sum(iS4_gt-iS4) < accuracy ;
res_partial(5) = sum(iS5_gt-iS5) < accuracy ;
res_partial(6) = sum(iS6_gt-iS6) < accuracy ;
res_partial(7) = sum(percentages1_gt-percentages1) < accuracy ;
res_partial(8) = sum(percentages2_gt-percentages2) < accuracy ;
res_partial(9) = sum(percentages3_gt-percentages3) < accuracy ;
res_partial(10) = sum(percentages4_gt-percentages4) < accuracy ;
res_partial(11) = sum(percentages5_gt-percentages5) < accuracy ;
res_partial(12) = sum(percentages6_gt-percentages6) < accuracy ;

res = all(res_partial);

%------------- END OF CODE --------------