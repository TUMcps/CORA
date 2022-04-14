function val = compIntersectionParam(W1,q1,W2,q2)
% compIntersectionParam - computes zero root of 'rootfnc' and returns the
% corresponding argument
%
% Syntax:  
%    val = compIntersectionParam(W1,q1,W2,q2)
%
% Inputs:
%    W1,q,W2,q2       - center and shape matrices of E1,E2 (see 'and')
%
% Outputs:
%    val - solution minimizing the volume of the parametrization in 'and'
%    (for details, see [1],[2])
%
% Example: 
%    -
% References: [1] Largely based on the Ellipsoidal Toolbox, see 
%                 https://de.mathworks.com/matlabcentral/fileexchange/21936-ellipsoidal-toolbox-et
%             [2] For detailed explanations, see:
%                 https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: and

% Author:       Victor Gassmann
% Written:      14-October-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
f = @(p) rootfnc(p,W1,q1,W2,q2);
val = fzero(f,0.5);
if (val < 0) || (val > 1) 
    if det(W1) > det(W2)
      val = 1;
    else
      val = 0;
    end
end
%------------- END OF CODE --------------