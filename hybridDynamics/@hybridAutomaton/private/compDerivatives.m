function compDerivatives(HA,options)
% compDerivatives - computes the derivatives of the flow equation for every
%    location with a nonlinear flow equation
%
% Syntax:
%    compDerivatives(HA,options)
%
% Inputs:
%    HA - hybridAutomaton object
%    options - reachability settings
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton/reach

% Authors:       Niklas Kochdumper
% Written:       20-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% loop over all locations
for i=1:length(HA.location)
   
    % read out location and corresponding flow equation
    loc = HA.location(i);
    sys = loc.contDynamics;
    
    % derivatives computation only required for nonlinear systems
    if isa(sys,'nonlinearSys') || isa(sys,'nonlinDASys') || ...
        isa(sys,'nonlinParamSys')
   
        % compute derivatives (generates files in models/auxiliary)
        derivatives(sys,options); 
    end
end

% ------------------------------ END OF CODE ------------------------------
