function [ output_args ] = test1( ~ )
% Example 2.2 FitzHugh-Nagumo system example from the article Non-linear Continuous Systems for
% Safety Verication(Benchmark Proposal) by Sogokon, Ghorbal, and Johnson
% Applied Verification for Continuous and Hybrid Systems
% /dot{x} = - x.^3 ./ 3 + x - y + 7/8
% /dot{y} = 2 / 25 .* ( x - 4 .* y ./ 5 + 7/10 )

x = interval(-1, -0.5)
y = interval(1, 1.5)
x1 = taylexp(x, 10);
y1 = taylexp(y, 10);

disp(' ')

tic
interval_x = - x.^3 ./ 3 + x - y + 7/8
t_interval_x = toc

tic
global_x = GInt2( x, y, '- x.^3 ./ 3 + x - y + 7/8', 0.01)
t_global_x = toc

tic
%tayl_roots_x = interval(- x1.^3 ./ 3 + x1 + 7/8) + interval( - y1)
tayl_int_x = interval(- x1.^3 ./ 3 + x1 + 7/8 - y1 )
t_tayl_int_x = toc

tic
tayl_gl_x = intervalGl2(- x1.^3 ./ 3 + x1 + 7/8 - y1, x, y, 0.001 )
t_tayl_gl_x = toc

disp(' ')

tic
interval_y = 2 / 25 .* ( x - 4 .* y ./ 5 + 7/10 )
t_interval_y = toc

tic
global_y = GInt2( x, y, '2 / 25 .* ( x - 4 .* y ./ 5 + 7/10 )', 0.01)
t_global_y = toc

tic
tayl_int_y = 2 / 25 .* ( interval( x1 - 4 .* y1 ./ 5 + 7/10 ))
t_tayl_int_y = toc

tic
tayl_gl_y = 2 / 25 .* ( intervalGl2( x1 - 4 .* y1 ./ 5 + 7/10 , x, y, 0.001))
t_tayl_gl_y = toc

end

