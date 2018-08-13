function completed = example_taylm_method1()
% updated: 29-June-2018, MA

syms a1 a2; % instantiate symbolic variables
s = [2 + 1.5*a1; 2.75 + 0.25*a2]; % create symbolic function
c = taylm(s, interval([-2;-3],[0;1]), 6)  % generate Taylor model

%example completed
completed = 1;

%------------- END OF CODE --------------