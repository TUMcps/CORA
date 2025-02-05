function [Obj]=stretch(Obj,factor)
% stretch - stretches the time constant T of the Markov chains 
% Pre:      Markov Chain object, factor
% Post:     Markov Chain object
%
% Syntax:
%    [Obj]=stretch(Obj,factor)
%
% Inputs:
%    ???
%
% Outputs:
%    ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       15-September-2006
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% ------------------------------ END OF CODE ------------------------------
