function [ output_args ] = test4( ~ )
% BRUSSELATOR or a chemical oscillator from CONTINUOUS-TIME CASE STUDIES (A)
% http://systems.cs.colorado.edu/research/cyberphysical/taylormodels/casestudies/
%
% /dot{x} = x.^2 .* y - 1.5 * x - x + 1
% /dot{y} = - x.^2 .* y + 1.5 * x

x = interval(0.9, 1)
y = interval(0, 0.1)
x1 = taylexp(x, 10);
y1 = taylexp(y, 10);

disp(' ')

tic
interval_x = x.^2 .* y - 1.5 * x - x + 1
t_interval_x = toc

tic
global_x = GInt2( x, y, 'x.^2 .* y - 1.5 * x - x + 1', 0.01)
t_global_x = toc

tic
tayl_int_x = interval( x1.^2 .* y1 - 1.5 * x1 - x1 + 1 )
t_tayl_int_x = toc

tic
tayl_gl_x = tayl2intGl2( x1.^2 .* y1 - 1.5 * x1 - x1 + 1, x, y, 0.01 )
t_tayl_gl_x = toc

disp(' ')

tic
interval_y = - x.^2 .* y + 1.5 * x
t_interval_y = toc

tic
global_y = GInt2( x, y, '- x.^2 .* y + 1.5 * x', 0.01)
t_global_y = toc

tic
tayl_int_y = interval( - x1.^2 .* y1 + 1.5 * x1 )
t_tayl_int_y = toc

tic
tayl_int_y = tayl2intGl2( - x1.^2 .* y1 + 1.5 * x1, x, y, 0.01 )
t_tayl_int_y = toc

end

