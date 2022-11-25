function [Obj]=stretch(Obj,factor)
% Purpose:  stretches the time constant T of the Markov chains 
% Pre:      Markov Chain object, factor
% Post:     Markov Chain object
% Tested:   15.09.06,MA


h = waitbar(0,'stretch T');

%T_T
T_pow{1}=Obj.T.T;
for i=2:factor
    T_pow{i}=T_pow{i-1}*Obj.T.T;
    %update waitbar
    waitbar(i/(2*factor),h);
end
Obj.T.T=T_pow{i};

%T_OT
T_temp=Obj.T.OT;
for i=2:factor
    T_temp=T_temp+Obj.T.OT*T_pow{i-1};
    %update waitbar
    waitbar((i+factor)/(2*factor),h);
end
Obj.T.OT=T_temp/factor;

%close waitbar
close(h);