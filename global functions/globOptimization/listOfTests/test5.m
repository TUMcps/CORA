function [ output_args ] = test5( ~ )
% LORENZ SYSTEM from CONTINUOUS-TIME CASE STUDIES (B)
% Lorenz System is a three variable dynamical system that can exhibit chaotic behavior
% http://systems.cs.colorado.edu/research/cyberphysical/taylormodels/casestudies/
%
% /dot{x} = 10 * (y - x)
% /dot{y} = x * (28 - z) - y
% /dot{z} = x * y - 8/3 * z

x = interval(14.8, 15.2)
y = interval(14.8, 15.2)
z = interval(35.8, 36.2)
x1 = taylexp(x, 10);
y1 = taylexp(y, 10);
z1 = taylexp(z, 10);

disp(' ')

tic
interval_x = 10 * (y - x)
t_interval_x = toc

tic
global_x = GInt3( x, y, z, '10 * (y - x)', 0.01)
t_global_x = toc

tic
tayl_int_x = interval( 10 * (y1 - x1) )
t_tayl_int_x = toc

tic
tayl_gl_x = tayl2intGl2( 10 * (y1 - x1), x, y, 0.01 )
t_tayl_gl_x = toc

disp(' ')

tic
interval_y = x .* (28 - z) - y
t_interval_y = toc

tic
global_y = GInt3( x, y, z, 'x .* (28 - z) - y', 0.01)
t_global_y = toc

tic
tayl_int_y = interval( x1 .* (28 - z1) - y1 )
t_tayl_int_y = toc

tic
tayl_int_y = tayl2intGl3( x1 .* (28 - z1) - y1, x, y, z, 0.01 )
t_tayl_int_y = toc

disp(' ')

tic
interval_z = x .* y - 8/3 * z
t_interval_z = toc

tic
global_z = GInt3( x, y, z, 'x .* y - 8/3 * z', 0.01)
t_global_z = toc

tic
tayl_int_z = interval( x1 .* y1 - 8/3 * z1 )
t_tayl_int_z = toc

tic
tayl_int_z = tayl2intGl3( x1 .* y1 - 8/3 * z1, x, y, z, 0.01 )
t_tayl_int_z = toc

end


