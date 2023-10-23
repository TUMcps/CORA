function Rint = observe_intersectionMethod_I(obj,R,y,options)
% observe_intersectionMethod_I - intersects the reachable set with
%    measurement strips according to intersection method I in [1]. 
%
% Syntax:
%    Rint = observe_intersectionMethod_I(obj,R,y,options)
%
% Inputs:
%    obj - discrete-time linear system object
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
% Written:       08-September-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init intersection zonotope
Rint = R;

for iStrip = 1:length(y)
    % in case a lambda value has been precomputed
    if isnumeric(options.intersectionTechnique)
        technique = options.intersectionTechnique(:,iStrip);
    else
        technique = options.intersectionTechnique;
    end
        
    % intersection of zonotope with strip
    Rint = intersectStrip(Rint,obj.C(iStrip,:),options.sigma(iStrip),...
        y(iStrip),technique);
end

% %% plot for debugging
% % convert strip to polytope
% C = [c{1}; -c{1}];
% d = [sigma{1} + y_strip{1}; sigma{1} - y_strip{1}];
% P = polytope(C,d);
% 
% % plot
% figure
% plot(R);
% hold on
% plot(P,[1 2],'k');
% plot(Rint,[1 2],'r--');

% ------------------------------ END OF CODE ------------------------------
