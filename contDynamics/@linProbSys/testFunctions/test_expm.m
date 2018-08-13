function [Exp_min,Exp_max,A_min,A_max]=test_expm(S,I,t)
% Purpose:  receive interval of matrix exponential by simulation 
% Pre:      single value matrix, interval matrix, time t
% Post:     simulated matrix exponential
% Tested:   19.01.07,MA

segments=[5,5,5,5];
nrOfSegments=prod(segments);
%init Exp_min, Exp_max, A
Exp_min=expm(S*t);
Exp_max=expm(S*t);
A=S;
A_min=S; A_max=S;
%------------------------
%get deltaI--------------
deltaI=I;
for iPos=1:length(segments)
    deltaI(iPos)=2*I(iPos)/(segments(iPos)-1);
end
%-----------------------
%compute system matrix A
for i=1:nrOfSegments
    seg=i2s(segments,i);
    for iPos=1:length(seg)
        A(iPos)=S(iPos)-I(iPos)+(seg(iPos)-1)*deltaI(iPos);
    end
    Exp=expm(A*t);
    for iPos=1:length(seg)
        %Exp_min/max
        if Exp(iPos)<Exp_min(iPos)
            Exp_min(iPos)=Exp(iPos);
        elseif Exp(iPos)>Exp_max(iPos)
            Exp_max(iPos)=Exp(iPos);
        end
        %A_min/max
        if A(iPos)<A_min(iPos)
            A_min(iPos)=A(iPos);
        elseif A(iPos)>A_max(iPos)
            A_max(iPos)=A(iPos);
        end
    end
end