function HA = spacecraft_levelSet()
% spacecraft_levelSet- spacecraft rendevous benchmark described in [1] with
%                      nonlinear guard sets
%
% Syntax:  
%    HA = spacecraft_levelSet()
%
% Inputs:
%    no
%
% Outputs:
%    HA - hybrid automaton object
%
% References: 
%   [1] N. Chan et al. "Verifying safety of an autonomous spacecraft 
%       rendezvous mission (Benchmark proposal)"  

% Author:       Niklas Kochdumper
% Written:      23-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Mode Approaching --------------------------------------------------------

%       x'  = vx 
%       y'  = vy 
%       vx' = (n^2 + K1_11/m_c)*x + (2*n + K1_14/m_c)*vy + K1_12/m_c * y + K1_13/m_c * vx + mu/r^2 - mu/sqrt((r+x)^2-y^2)^3 * (r+x)
%       vy' = (n^2 + K1_22/m_c)*y + (K1_23/m_c -2*n)*vx + K1_21/m_c * x + K1_24/m_c * vy - mu/sqrt((r+x)^2-y^2)^3 * y

% system dynamics
dynamics = nonlinearSys(@dynamics_approaching); 

% invariant:  x^2 + y^2 >= 100^2
syms x y vx vy;
inv = levelSet(-x^2 - y^2 + 100^2,[x;y;vx;vy],'<='); 

% transition: x^2 + y^2 == 100^2
resetA = eye(4);
resetc = zeros(4,1);
reset = struct('A', resetA, 'c', resetc);

guard = levelSet(x^2 + y^2 - 100^2,[x;y;vx;vy],'==');

trans = transition(guard, reset, 2);

% location
loc(1) = location('S1', inv, trans, dynamics);


% Mode Rendezvous Attempt -------------------------------------------------

%       x'  = vx
%       y'  = vy 
%       vx' = (n^2 + K2_11/m_c)*x + (2*n + K2_14/m_c)*vy + K2_12/m_c * y + K2_13/m_c * vx + mu/r^2 - mu/sqrt((r+x)^2-y^2)^3 * (r+x)
%       vy' = (n^2 + K2_22/m_c)*y + (K2_23/m_c -2*n)*vx + K2_21/m_c * x + K2_24/m_c * vy - mu/sqrt((r+x)^2-y^2)^3 * y
       
% system dynamics
dynamics = nonlinearSys(@dynamics_attempt); 

% invariant: x^2 + y^2 < 100^2
syms x y vx vy;
inv = levelSet(x^2 + y^2 - 100^2,[x;y;vx;vy],'<'); 

% location
trans = transition();
loc(2) = location('S2', inv, trans, dynamics);



% Hybrid Automaton --------------------------------------------------------

HA = hybridAutomaton(loc);


end

%------------- END OF CODE --------------