function [ output_args ] = test6( ~ )
% Hybrid case studies in the first mode l_1 from Hybrid
% Lorenz System is a three variable dynamical system that can exhibit chaotic behavior
% http://systems.cs.colorado.edu/research/cyberphysical/taylormodels/casestudies/
%
% /dot{x} = -9*(x - 2) - 7 * (y + 2) + (z - 1) + 0.2*(x - 2)*(y + 2) + 0.1 * (y + 2) .* (z - 1) + 0.1 * (x - 2) .* (z - 1) + 0.5 * (z - 1).^2 
% /dot{y} = 6*(x - 2) + 4*(y + 2) + (z - 1)
% /dot{z} = 3 *(x - 2) + 2 *(y + 2) - 2.5*(z - 1)

x = interval(3, 3.5)
y = interval(-3, -2.5)
z = interval(0.7, 1.3)
x1 = taylexp(x, 10);
y1 = taylexp(y, 10);
z1 = taylexp(z, 10);

disp(' ')

tic
interval_x = -9*(x - 2) - 7 * (y + 2) + (z - 1) + 0.2 * (x - 2).*(y + 2) + 0.1 * (y + 2) .* (z - 1) + 0.1 * (x - 2) .* (z - 1) + 0.5 * (z - 1).^2
t_interval_x = toc

tic
global_x = GInt3( x, y, z, '-9*(x - 2) - 7 * (y + 2) + (z - 1) + 0.2 * (x - 2).*(y + 2) + 0.1 * (y + 2) .* (z - 1) + 0.1 * (x - 2) .* (z - 1) + 0.5 * (z - 1).^2', 0.01)
t_global_x = toc

tic
tayl_int_x = interval( -9*(x1 - 2) - 7 * (y1 + 2) + (z1 - 1) + 0.2 * (x1 - 2).*(y1 + 2) + 0.1 * (y1 + 2) .* (z1 - 1) + 0.1 * (x1 - 2) .* (z1 - 1) + 0.5 * (z1 - 1).^2 ) 
t_tayl_int_x = toc

tic
tayl_gl_x = tayl2intGl3( -9*(x1 - 2) - 7 * (y1 + 2) + (z1 - 1) + 0.2 * (x1 - 2).*(y1 + 2) + 0.1 * (y1 + 2) .* (z1 - 1) + 0.1 * (x1 - 2) .* (z1 - 1) + 0.5 * (z1 - 1).^2, x, y, z, 0.01 ) 
t_tayl_gl_x = toc

disp(' ')

tic
interval_y = 6*(x - 2) + 4*(y + 2) + (z - 1)
t_interval_y = toc

tic
global_y = GInt3( x, y, z, '6*(x - 2) + 4*(y + 2) + (z - 1)', 0.01)
t_global_y = toc

tic
tayl_int_y = interval( 6*(x1 - 2) + 4*(y1 + 2) + (z1 - 1) )
t_tayl_int_y = toc

tic
tayl_gl_y = tayl2intGl3( 6*(x1 - 2) + 4*(y1 + 2) + (z1 - 1), x, y, z, 0.01 )
t_tayl_gl_y = toc

disp(' ')

tic
interval_z = 3 *(x - 2) + 2 *(y + 2) - 2.5*(z - 1)
t_interval_z = toc

tic
global_z = GInt3( x, y, z, '3 *(x - 2) + 2 *(y + 2) - 2.5*(z - 1)', 0.01)
t_global_z = toc

tic
tayl_int_z = interval( 3 *(x1 - 2) + 2 *(y1 + 2) - 2.5*(z1 - 1) )
t_tayl_int_z = toc

tic
tayl_gl_z = tayl2intGl3( 3 *(x1 - 2) + 2 *(y1 + 2) - 2.5*(z1 - 1), x, y, z, 0.01 )
t_tayl_gl_z = toc

end

