function [ output_args ] = test8( ~ )
% INSULIN INFUSION CONTROL IN DIABETIC PATIENTS from Hybrid
% dynamics of insulin glucose in a type I diabetic patient is modeled by 
% the Bergman minimal model with three state variables (G,I,X) wherein G 
% is plasma glucose concentration above the basal value G_B, I is 
% the plasma insulin concentration above the basal value I_B, and X is 
% the insulin concentration in an interstitial chamber.
% http://systems.cs.colorado.edu/research/cyberphysical/taylormodels/casestudies/
%
% /dot{x} = -0.01 * y - x .* (y + 4.5)
% /dot{y} = -0.025 * x + 0.000013 * z
% /dot{z} = -0.093 * (z + 15) + 0.083 * 8.333

x = interval(-0.1, 0.1)
y = interval(-2, 2)
z = interval(-0.1, 0.1)
x1 = taylexp(x, 10);
y1 = taylexp(y, 10);
z1 = taylexp(z, 10);

disp(' ')

tic
interval_x = -0.01 * y - x .* (y + 4.5)
t_interval_x = toc

tic
global_x = GInt2( x, y, '-0.01 * y - x .* (y + 4.5)', 0.01)
t_global_x = toc

tic
tayl_int_x = interval( -0.01 * y1 - x1 .* (y1 + 4.5) )
t_tayl_int_x = toc

tic
tayl_gl_x = tayl2intGl2( -0.01 * y1 - x1 .* (y1 + 4.5), x, y, 0.01 )
t_tayl_gl_x = toc

disp(' ')

tic
interval_y = -0.025 * x + 0.000013 * z
t_interval_y = toc

tic
global_y = GInt2( x, z, '-0.025 * x + 0.000013 * y', 0.01)
t_global_y = toc

tic
tayl_int_y = interval( -0.025 * x1 + 0.000013 * y1 )
t_tayl_int_y = toc

tic
tayl_gl_y = tayl2intGl2( -0.025 * x1 + 0.000013 * y1, x, y, 0.01 )
t_tayl_gl_y = toc

disp(' ')

tic
interval_z = -0.093 * (z + 15) + 0.083 * 8.333
t_interval_z = toc

tic
global_z = GInt1( z, '-0.093 * (x + 15) + 0.083 * 8.333', 0.01)
t_global_z = toc

tic
tayl_int_z = interval( -0.093 * (x1 + 15) + 0.083 * 8.333 )
t_tayl_int_z = toc

tic
tayl_int_z = tayl2intGl1( -0.093 * (x1 + 15) + 0.083 * 8.333, z, 0.01 )
t_tayl_int_z = toc

end

