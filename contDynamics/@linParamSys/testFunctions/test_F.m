function [F_min,F_max,A_min,A_max]=test_F(S,I,r)
% Purpose:  receive interval of matrix exponential by simulation 
% Pre:      single value matrix, interval matrix, time t
% Post:     simulated matrix exponential
% Tested:   19.01.07,MA

times=10;
delta_t=r/times;
segments=[5,5,5,5];
nrOfSegments=prod(segments);
dim=length(S);
%init Exp_min, Exp_max, A
A=S;
A_min=S; A_max=S;
%------------------------
%get deltaI--------------
deltaI=I;
for iPos=1:length(segments)
    deltaI(iPos)=2*I(iPos)/(segments(iPos)-1);
end
%-----------------------
%compute times
for iTime=1:times
    t=delta_t*iTime;
    %compute system matrix A
    for i=1:nrOfSegments
        seg=i2s(segments,i);
        for iPos=1:length(seg)
            A(iPos)=S(iPos)-I(iPos)+(seg(iPos)-1)*deltaI(iPos);
        end
        Exp_t=expm(A*t);
        Exp_r=expm(A*r);
        F=Exp_t-eye(dim)-t/r*(Exp_r-eye(dim));
        
        if (iTime==1) & (i==1)
            F_min=F;
            F_max=F;
        else
            for iPos=1:length(seg)
                %Exp_min/max
                if F(iPos)<F_min(iPos)
                    F_min(iPos)=F(iPos);
                elseif F(iPos)>F_max(iPos)
                    F_max(iPos)=F(iPos);
                end
                %A_min/max
                if A(iPos)<A_min(iPos)
                    A_min(iPos)=A(iPos);
                elseif A(iPos)>A_max(iPos)
                    A_max(iPos)=A(iPos);
                end
            end
        end
    end
end