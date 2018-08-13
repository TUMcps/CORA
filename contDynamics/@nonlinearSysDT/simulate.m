function [obj,k,x] = simulate(obj,options,kstart,kfinal,x0)
% simulate - simulates the system within a location
%
% Syntax:  
%    [t,x,index] = simulate(obj,kstart,kfinal,x0,options)
%
% Inputs:
%    obj - nonlinearSysDT object
%    kstart - start index
%    kfinal - final index
%    x0 - initial state 
%    options - structure containing the algorithm options
%
% Outputs:
%    obj - nonlinearSysDT object
%    k - index vector
%    x - state vector
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      22-August-2012
% Last update:  29-January-2018 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

%set initial state
x(:,1) = x0;
k(1) = kstart;

%loop
for i = kstart:kfinal
    
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,iSet);
    end
    
    x(:,i+1) = obj.mFile(0,x(:,i),options.uTrans,options.timeStep);
    
    k(i+1) = i;
end
    
%------------- END OF CODE --------------