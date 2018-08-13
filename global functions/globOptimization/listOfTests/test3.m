function [ output_args ] = test3( ~ )
% Example 2.3 from the article Non-linear Continuous Systems for
% Safety Verication(Benchmark Proposal) by Sogokon, Ghorbal, and Johnson
% Applied Verification for Continuous and Hybrid Systems
%
% Namely the ODE from Freddy Dumortier, Jaume Llibre, and Joan C. Artes.
% Qualitative  Theory  of  Planar  Diferential Systems.  Springer, 2006
%
% /dot{x1} = -42 * x1^7 + 68 * x1^6 * y1 - 46 * x1^5 .* y1 + 258 * x1.^4 .* y1 + 156 * x1.^3 .* y1 + 50 * x1.^2 .* y1 + 20 * x1 .* y.^6 - 8 * y1^7
% /dot{y1} = y1 .* ( 1110 .* x1^6 - 220 .* x1^5 .* y1 - 3182 .* x1^4 .* y1 + 478 .* x1^3 .* y1^3 + 487 .* x1^2 .* y1^4 - 102 .* x1 .* y1^5 - 12 .* y1^6 )

x = interval(-1, -3/4)
y = interval(1, 1.5)
x1 = taylexp(x, 10);
y1 = taylexp(y, 10);

disp(' ')

tic
interval_x = -42 * x.^7 + 68 * x.^6 .* y - 46 .* x.^5 .* y + 258 * x.^4 .* y + 156 * x.^3 .* y + 50 * x.^2 .* y + 20 * x .* y.^6 - 8 * y.^7
t_interval_x = toc

tic
global_x = GInt2( x, y, '-42 * x.^7 + 68 * x.^6 .* y - 46 .* x.^5 .* y + 258 * x.^4 .* y + 156 * x.^3 .* y + 50 * x.^2 .* y + 20 * x .* y.^6 - 8 * y.^7', 0.01)
t_global_x = toc

tic
tayl_int_x = interval(-42 * x1.^7 + 68 * x1.^6 .* y1 - 46 .* x1.^5 .* y1 + 258 * x1.^4 .* y1 + 156 * x1.^3 .* y1 + 50 * x1.^2 .* y1 + 20 * x1 .* y1.^6 - 8 * y1.^7)
t_tayl_int_x = toc

tic
tayl_gl_x = tayl2intGl2(-42 * x1.^7 + 68 * x1.^6 .* y1 - 46 .* x1.^5 .* y1 + 258 * x1.^4 .* y1 + 156 * x1.^3 .* y1 + 50 * x1.^2 .* y1 + 20 * x1 .* y1.^6 - 8 * y1.^7, x, y, 0.01)
t_tayl_gl_x = toc

disp(' ')

tic
interval_y = y .* ( 1110 .* x.^6 - 220 .* x.^5 .* y - 3182 .* x.^4 .* y + 478 .* x.^3 .* y.^3 + 487 .* x.^2 .* y.^4 - 102 .* x .* y.^5 - 12 .* y.^6 )
t_interval_y = toc

tic
global_y = GInt2( x, y, 'y .* ( 1110 .* x.^6 - 220 .* x.^5 .* y - 3182 .* x.^4 .* y + 478 .* x.^3 .* y.^3 + 487 .* x.^2 .* y.^4 - 102 .* x .* y.^5 - 12 .* y.^6 )', 0.01)
t_global_y = toc

tic
tayl_int_y = interval( y1 .* ( 1110 .* x1.^6 - 220 .* x1.^5 .* y1 - 3182 .* x1.^4 .* y1 + 478 .* x1.^3 .* y1.^3 + 487 .* x1.^2 .* y1.^4 - 102 .* x1 .* y1.^5 - 12 .* y1.^6 ) )
t_tayl_int_y = toc

tic
tayl_gl_y = tayl2intGl2( y1 .* ( 1110 .* x1.^6 - 220 .* x1.^5 .* y1 - 3182 .* x1.^4 .* y1 + 478 .* x1.^3 .* y1.^3 + 487 .* x1.^2 .* y1.^4 - 102 .* x1 .* y1.^5 - 12 .* y1.^6 ), x, y, 0.01 )
t_tayl_gl_y = toc

end

