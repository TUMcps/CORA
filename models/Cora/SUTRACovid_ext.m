function f = SUTRACovid_ext(x,u)
% SUTRACovid_ext - SUTRA Covid disease spread benchmark [1] with
%    extended state vector (parameters uncertain but constant)
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
% parameters: beta, gamma -> extended states

f(1,1) = -x(7)*x(1) * (x(3)+x(4));
f(2,1) = -x(7)*x(2) * (x(3)+x(4));
f(3,1) = x(7)*x(1) * (x(3)+x(4)) - x(8)*x(3);
f(4,1) = x(7)*x(2) * (x(3)+x(4)) - x(8)*x(4);
f(5,1) = x(8)*x(3);
f(6,1) = x(8)*x(4);
f(7,1) = 0;
f(8,1) = 0;

end

%------------- END OF CODE --------------
