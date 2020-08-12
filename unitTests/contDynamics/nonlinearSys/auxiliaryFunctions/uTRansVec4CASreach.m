function uTransVec = uTRansVec4CASreach()
% updated: 05-August-2013, MA


%load data
load('linearized_controller_09_double_lane_change_jy(0)=0,jy(1)=0.mat');

%write values for each point in time
for i = 1:length(R)
    %R
    uTransVec(1,i) = R{i}(1,1);
    uTransVec(2,i) = R{i}(1,2);
    uTransVec(3,i) = R{i}(1,3);
    uTransVec(4,i) = R{i}(1,4);
    uTransVec(5,i) = R{i}(1,5);
    uTransVec(6,i) = R{i}(1,6);
    uTransVec(7,i) = R{i}(1,7);
    uTransVec(8,i) = R{i}(1,8);
    uTransVec(9,i) = R{i}(2,1);
    uTransVec(10,i) = R{i}(2,2);
    uTransVec(11,i) = R{i}(2,3);
    uTransVec(12,i) = R{i}(2,4);
    uTransVec(13,i) = R{i}(2,5);
    uTransVec(14,i) = R{i}(2,6);
    uTransVec(15,i) = R{i}(2,7);
    uTransVec(16,i) = R{i}(2,8);
    
    %Xn
    uTransVec(17,i) = Xn{i}(1);
    uTransVec(18,i) = Xn{i}(2);
    uTransVec(19,i) = Xn{i}(3);
    uTransVec(20,i) = Xn{i}(4);
    uTransVec(21,i) = Xn{i}(5);
    uTransVec(22,i) = Xn{i}(6);
    uTransVec(23,i) = Xn{i}(7);
    uTransVec(24,i) = Xn{i}(8);
    
    %W
    uTransVec(25,i) = W{i}(1);
    uTransVec(26,i) = W{i}(2);
end



