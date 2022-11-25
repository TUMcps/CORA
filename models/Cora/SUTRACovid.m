function f = SUTRACovid(x,u,p)
% SUTRACovid - SUTRA Covid disease spread benchmark [1]
%
% Syntax:  
%    f = SUTRACovid(x,u,p)
%
% Inputs:
%    x - state vector
%    u - input vector
%    p - parameter vector
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

f(1,1) = -p(1)*x(1) * (x(3)+x(4));
f(2,1) = -p(1)*x(2) * (x(3)+x(4));
f(3,1) = p(1)*x(1) * (x(3)+x(4)) - p(2)*x(3);
f(4,1) = p(1)*x(2) * (x(3)+x(4)) - p(2)*x(4);
f(5,1) = p(2)*x(3);
f(6,1) = p(2)*x(4);

end

%------------- END OF CODE --------------
