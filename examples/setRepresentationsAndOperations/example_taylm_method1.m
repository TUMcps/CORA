function completed = example_taylm_method1()
% updated: 29-June-2018, MA

a1 = interval(-1, 2); % generate a scalar interval [-1,2]
b1 = taylm(a1, 6);	  % generate a scalar Taylor model of order 6
a2 = interval(2, 3);  % generate a scalar interval [2,3]
b2 = taylm(a2, 6);	  % generate a scalar Taylor model of order 6
c = [b1; b2]		  % generate a row of Taylor models

%example completed
completed = 1;

%------------- END OF CODE --------------