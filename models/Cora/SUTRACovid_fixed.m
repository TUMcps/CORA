function f = SUTRACovid_fixed(x,u)
% SUTRACovid_fixed - SUTRA Covid disease spread benchmark [1] without parameters
%
% Syntax:  
%    f = SUTRACovid_fixed(x,u,p)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - time derivate of the system state
%
% Reference:
%   [1] National Supermodel Committee. Indian supermodel for covid-19 pandemic.
%       https://www.iith.ac.in/~m_vidyasagar/arXiv/Super-Model.pdf, May 2021.

% Author:       Mark Wetzlinger
% Written:      18-May-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% states:     S_A, S_I, A, I, R_I, R_A
% parameters: beta = 0.25, gamma = 0.02

f(1,1) = -0.25*x(1) * (x(3)+x(4));
f(2,1) = -0.25*x(2) * (x(3)+x(4));
f(3,1) = 0.25*x(1) * (x(3)+x(4)) - 0.02*x(3);
f(4,1) = 0.25*x(2) * (x(3)+x(4)) - 0.02*x(4);
f(5,1) = 0.02*x(3);
f(6,1) = 0.02*x(4);

end

%------------- END OF CODE --------------
