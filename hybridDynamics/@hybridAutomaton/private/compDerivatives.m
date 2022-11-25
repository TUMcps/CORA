function compDerivatives(obj,options)
% compDerivatives - compute the derivatives of the dynamic function for
%                   nonlinear systems
%
% Syntax:  
%    compDerivatives(obj,options)
%
% Inputs:
%    obj - hybrid automaton object
%    options - options for the computation of the reachable set
%
% Outputs:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton/reach

% Author:       Niklas Kochdumper
% Written:      20-May-2020
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    for i = 1:length(obj.location)
       
        loc = obj.location{i};
        sys = loc.contDynamics;
        
        if isa(sys,'nonlinearSys') || isa(sys,'nonlinDASys') || ...
           isa(sys,'nonlinParamSys')
       
           derivatives(sys,options); 
        end
    end

%------------- END OF CODE --------------