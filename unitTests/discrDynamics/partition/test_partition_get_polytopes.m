function res = test_partition_get_polytopes()
% test_partition_get_polytopes - unit test cell polytopes
%
% Syntax:  
%    res = test_partition_get_polytopes()
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


%setup partitions
threeDimField=partition([0,10; -3,3; 0,1],[5;10;3]);
twoDimField=partition([0,10; -3,3],[5;10]);

% check that cellPolytopes works, 3DOF
Pols1 = cellPolytopes(twoDimField,1:nrOfCells(twoDimField));
Pols2 = cellPolytopes(twoDimField);

if length(Pols1)==length(Pols2)
    res1 = (norm(volume(Pols1{3}) - volume(Pols2{3}))<1e-9);
else
    res1 = 0;
end


% check that cellPolytopes works, 3DOF
Pols1 = cellPolytopes(threeDimField,1:nrOfCells(threeDimField));
Pols2 = cellPolytopes(threeDimField);

if length(Pols1)==length(Pols2)
    res2 = (norm(volume(Pols1{3}) - volume(Pols2{3}))<1e-9);
else
    res2 = 0;
end


res = res1&&res2;
% 
% cellPolytopes(threeDimField,[1 5 3])
% cellPolytopes(threeDimField)
% cellPolytopes(threeDimField,[1 5 3])
% cellPolytopes(threeDimField)
% P = mptPolytope([2 0 0.3;4 2 0.6;1 1 0.5; 1 1 0.1]);
% intersectingSegments(threeDimField,P)
% [iS,percentages] = exactIntersectingCells(threeDimField,P)
% plot(threeDimField,exactIntersectingCells(threeDimField,P))
% hold on
% plot(P)
% 
% %partition with 

%------------- END OF CODE --------------