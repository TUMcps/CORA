function [obj]=reach_timeStep(obj,Zinit,T,order)
% reach - computes the reachable set of the specified linear interval
% system for a initial zonotope set Z_init and a certain time horizon T
%
% Syntax:  
%    [obj]=reach(obj,Zinit,T)
%
% Inputs:
%    obj - linear interval system object
%    Zinit - initial zonotope set
%    T - time horizon for the computation of reachable sets
%    order - specifies maximum order of zonotopes
%
% Outputs:
%    obj - linear interval system object
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author:       Matthias Althoff
% Written:      23-January-2007 
% Last update:  30-April-2007
%               15-June-2016
%               25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

%------------- BEGIN CODE --------------

%load data from object structure
eAt=obj.taylor.eAt;
F=obj.taylor.F;
inputSol=obj.taylor.input;
accuracy=obj.accuracy;
r=obj.r;
dim=obj.dim;

%first time step
R=eAt*Zinit;
%work around for enclose:
%Zinit=Zinit+0*zonotope(zeros(dim,dim+1));
%convex hull
%R=enclose(Zinit,R);
%convert inputSol to a zonotope called Zu
Zu=zonotope(interval(infimum(inputSol),supremum(inputSol)));
%enclarge due to time interval and input
%R=R+F*Zinit+Zu;
R=R+Zu;
%reduce zonotope
Z{1}=reduce(R,'girard',order);

%initialize time and counter i
time=0; i=0;

%time horizon not reached
while time<T
    %increase counter, time
    i=i+1; time=time+r;
    %next zonotope element (time point, no input)
    R=eAt*Z{i};
    %bloating due to input
    R=R+Zu;
    %reduce
    Z{i+1}=reduce(R,'girard',order);
end
%write zonotopes to object structure
obj.reachSet=Z;
    
%------------- END OF CODE --------------