function [ output_args ] = test7( ~ )
% Hybrid studies case in the second mode l_2 from Hybrid
% Lorenz System is a three variable dynamical system that can exhibit chaotic behavior
% http://systems.cs.colorado.edu/research/cyberphysical/taylormodels/casestudies/
%
% /dot{x} = 2.2 * x + 3.6 * y + 3.9 * z
% /dot{y} = 3 * x + 2.4 * y + 3.4 * z - 0.01 * x.^2
% /dot{z} = -5 * x - 5.4 * y - 6.7 * z

x = interval(3, 3.5)
y = interval(-3, -2.5)
z = interval(0.7, 1.3)
x1 = taylexp(x, 10);
y1 = taylexp(y, 10);
z1 = taylexp(z, 10);

disp(' ')

tic
interval_x = 2.2 * x + 3.6 * y + 3.9 * z
t_interval_x = toc

tic
global_x = GInt3( x, y, z, '2.2 * x + 3.6 * y + 3.9 * z', 0.01)
t_global_x = toc

tic
tayl_int_x = interval( 2.2 * x1 + 3.6 * y1 + 3.9 * z1 )
t_tayl_int_x = toc

tic
tayl_gl_x = tayl2intGl3( 2.2 * x1 + 3.6 * y1 + 3.9 * z1, x, y, z, 0.01 )
t_tayl_gl_x = toc

disp(' ')

tic
interval_y = 3 * x + 2.4 * y + 3.4 * z - 0.01 * x.^2
t_interval_y = toc

tic
global_y = GInt3( x, y, z, '3 * x + 2.4 * y + 3.4 * z - 0.01 * x.^2', 0.01)
t_global_y = toc

tic
tayl_int_y = interval( 3 * x1 + 2.4 * y1 + 3.4 * z1 - 0.01 * x1.^2 )
t_tayl_int_y = toc

tic
tayl_gl_y = tayl2intGl3( 3 * x1 + 2.4 * y1 + 3.4 * z1 - 0.01 * x1.^2, x, y, z, 0.01 )
t_tayl_gl_y = toc

disp(' ')

tic
interval_z = -5 * x - 5.4 * y - 6.7 * z
t_interval_z = toc

tic
global_z = GInt3( x, y, z, '-5 * x - 5.4 * y - 6.7 * z', 0.01)
t_global_z = toc

tic
tayl_int_z = interval( -5 * x1 - 5.4 * y1 - 6.7 * z1 )
t_tayl_int_z = toc

tic
tayl_int_z = tayl2intGl3( -5 * x1 - 5.4 * y1 - 6.7 * z1, x, y, z, 0.01 )
t_tayl_int_z = toc

end

