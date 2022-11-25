function f = SUTRACovid_inputs(x,u)
% SUTRACovid_inputs - SUTRA Covid disease spread benchmark [1], where
%    the uncertain parameters are modelled as inputs
%
% Syntax:  
%    f = SUTRACovid(x,u,p)
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
% Written:      17-May-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% states:     S_A, S_I, A, I, R_I, R_A
% parameters: beta, gamma

f(1,1) = -u(1)*x(1) * (x(3)+x(4));
f(2,1) = -u(1)*x(2) * (x(3)+x(4));
f(3,1) = u(1)*x(1) * (x(3)+x(4)) - u(2)*x(3);
f(4,1) = u(1)*x(2) * (x(3)+x(4)) - u(2)*x(4);
f(5,1) = u(2)*x(3);
f(6,1) = u(2)*x(4);

end

%------------- END OF CODE --------------
