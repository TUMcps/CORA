function completed = example_taylm_method2()
% updated: 29-June-2018, MA

a = interval([-1;2], [2;3]);	% generate an interval vector [[-1,2]; [2,3]]
c = taylm(a, 6, {'a1';'a2'}) 	% generate a Taylor model (order 6) with variable names a1 and a2

%example completed
completed = 1;

%------------- END OF CODE --------------