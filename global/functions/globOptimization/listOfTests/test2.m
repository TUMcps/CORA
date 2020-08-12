function [ output_args ] = test2( ~ )
% Example 2.1 Non-linear example from the article Non-linear Continuous Systems for
% Safety Verication(Benchmark Proposal) by Sogokon, Ghorbal, and Johnson
% Applied Verification for Continuous and Hybrid Systems
% /dot{x} = 2 * x - x .* y
% /dot{y} = 2 * x .^ 2 - y

x = interval(-1, -0.5)
y = interval(1, 1.5)
x1 = taylexp(x, 10);
y1 = taylexp(y, 10);

disp(' ')

tic
interval_x = 2 * x - x .* y
t_interval_x = toc

tic
global_x = GInt2( x, y, '2 * x - x .* y', 0.01)
t_global_x = toc

tic
tayl_int_x = interval(2 * x1 - x1 .* y1 )
t_tayl_int_x = toc

tic
tayl_gl_x = intervalGl2(2 * x1 - x1 .* y1, x, y, 0.001 )
t_tayl_gl_x = toc

disp(' ')

tic
interval_y = 2 * x .^ 2 - y
t_interval_y = toc

tic
global_y = GInt2( x, y, '2 * x .^ 2 - y', 0.01)
t_global_y = toc

tic
tayl_int_y = interval( 2 * x1 ^ 2  - y1 ) 
t_tayl_int_y = toc

tic
tayl_gl_y = intervalGl2( 2 * x1 ^ 2  - y1, x, y, 0.001 ) 
t_tayl_gl_y = toc

end

