function vol = modVolume(P)
% modVolume - reliable computation of the volume of a polytope; if volume
%    cannot be computed, it is overapproximated by an interval hull
%
% Syntax:
%    vol = modVolume(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    vol - exact or overapproximated volume of a polytope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       29-October-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    vol=volume(P);
catch
    %overapproximate polytope by an interval hull
    V=extreme(P)';
    IH=interval.enclosePoints(V);
    vol=volume(IH);        
    disp('Volume had to be overapproximated')
end

% ------------------------------ END OF CODE ------------------------------
