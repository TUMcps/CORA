function Rint = priv_observe_intersectionMethod_II(sys,R,y,options)
% priv_observe_intersectionMethod_II - intersects the reachable set with
%    measurement strips according to intersection method II in [1]. 
%
% Syntax:
%    Rint = priv_observe_intersectionMethod_II(sys,R,y,options)
%
% Inputs:
%    sys - discrete-time linear system object
%    R - reachable set
%    y - current measurement
%    options - options for the computation
%
% Outputs:
%    Rint - resulting zonotope after intersections with strips
%
% Example: 
%    -
%
% Reference:
%    [1] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
%        for Guaranteed State Estimation of Linear Disturbed Systems, 
%        in preparation.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   02-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% refresh system matrix if system is nonlinear and the intersection method
% is 'wang-FRad'
if isa(sys,'nonlinearSysDT') && strcmp(options.intersectionTechnique.method,'wang-FRad')
    options.intersectionTechnique.A = sys.jacobian(center(R),[]); % system matrix
end

% intersection of zonotope with strips
Rint = intersectStrip(R,sys.C,options.sigma,y,options.intersectionTechnique);

% ------------------------------ END OF CODE ------------------------------
