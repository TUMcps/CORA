function completed = example_converter_powerSystem2cora_IEEE30()
% example_converter_powerSystem2cora_IEEE30 - example for creating the IEEE
%    30 bus power system benchmark, can be found in [1, Sec. VII].
%
% Syntax:
%    completed = example_converter_powerSystem2cora_IEEE30()
%
% Inputs:
%    true
%
% Outputs:
%    completed - true/false 
%
% References:
%    [1] M. Althoff, "Formal and Compositional Analysis of Power Systems 
%        using Reachable Sets", IEEE Transactions on Power Systems 29 (5), 
%        2014, 2270-2280

% Authors:       Matthias Althoff
% Written:       14-April-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% create CORA models
% full system
powerSystem2cora('IEEE30')
% subsystem 1
powerSystem2cora('IEEE30_sub1')
% subsystem 2
powerSystem2cora('IEEE30_sub2')

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
