function res = testLongDuration_mptPolytope_in
% testLongDuration_mptPolytope_in - unit test function for containment
%    check
%
% Syntax:  
%    res = testLongDuration_mptPolytope_in()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      26-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% Random Test -------------------------------------------------------------

for i = 2:2:4
   
    % create random polytope
    P = mptPolytope.generateRandom(i);
    
    % check if randomly generated points are inside
    Y = randPoint(P,2*i);
    if ~in(P,Y,1e-12)
        res = false;
        break;
    end
    
    % generate zonotope
    Z = zonotope.generateRandom(i);
    Ei = ellipsoid(Z,'i:norm');
    Pz = polytope(Z);
    % check if Ei is in Pz
    if ~in(Pz,Ei,Ei.TOL)
        res = false;
        break;
    end
end

if res
    disp('testLongDuration_mptPolytope_in successful');
else
    disp('testLongDuration_mptPolytope_in failed');
end
    