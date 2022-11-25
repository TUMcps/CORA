function dx = vehicle(x, u)
    % all motion states
    v_x = x(1);
    v_y = x(2);
    d_psi = x(3);
    delta_1 = x(4);
    delta_2 = x(5);
    delta_3 = x(6);
    delta_4 = x(7);
    omega_1 = x(8);
    omega_2 = x(9);
    omega_3 = x(10);
    omega_4 = x(11);
    % velocity of each wheel center
    G_trans = [1, 0, -0.7;
                   0, 1,  0.3;
                   1, 0,  0.3;
                   0, 1,  0.7;
                   1, 0, -0.3;
                   0, 1, -0.3;
                   1, 0,  0.7;
                   0, 1, -0.3];
    v_xy = G_trans * [v_x; v_y; d_psi];
    % wheel slip
    s1_x = (0.3 * omega_1 * cos(delta_1) - v_xy(1)) / sqrt(v_xy(1)^2+v_xy(2)^2);
    s1_y = (0.3 * omega_1 * sin(delta_1) - v_xy(2)) / sqrt(v_xy(1)^2+v_xy(2)^2);
    s2_x = (0.3 * omega_2 * cos(delta_2) - v_xy(3)) / sqrt(v_xy(3)^2+v_xy(4)^2);
    s2_y = (0.3 * omega_2 * sin(delta_2) - v_xy(4)) / sqrt(v_xy(3)^2+v_xy(4)^2);
    s3_x = (0.3 * omega_3 * cos(delta_3) - v_xy(5)) / sqrt(v_xy(5)^2+v_xy(6)^2);
    s3_y = (0.3 * omega_3 * sin(delta_3) - v_xy(6)) / sqrt(v_xy(5)^2+v_xy(6)^2);
    s4_x = (0.3 * omega_4 * cos(delta_4) - v_xy(7)) / sqrt(v_xy(7)^2+v_xy(8)^2);
    s4_y = (0.3 * omega_4 * sin(delta_4) - v_xy(8)) / sqrt(v_xy(7)^2+v_xy(8)^2);
    % simplified effective force of wheel
    Fx_1 = 800 * s1_x;
    Fy_1 = 800 * s1_y;
    Fx_2 = 800 * s2_x;
    Fy_2 = 800 * s2_y;
    Fx_3 = 800 * s3_x;
    Fy_3 = 800 * s3_y;
    Fx_4 = 800 * s4_x;
    Fy_4 = 800 * s4_y;
    % effective force of the vehicle center
    F_H = transpose(G_trans) * [Fx_1; Fy_1; Fx_2; Fy_2; Fx_3; Fy_3; Fx_4; Fy_4];
    % dyamics of the vehicle center
    a_H = [1/500, 0, 0; 0, 1/500, 0; 0, 0, 1/300] * F_H;
    % velocity of the vehicle center and nonlinear model
    dx(1, 1) = v_y * d_psi + a_H(1);
    dx(2, 1) = -v_x * d_psi + a_H(2);
    dx(3, 1) = a_H(3);
    dx(4, 1) = 0;
    dx(5, 1) = 0;
    dx(6, 1) = 0;
    dx(7, 1) = 0;
    dx(8, 1) = 0;
    dx(9, 1) = 0;
    dx(10, 1) = 0;
    dx(11, 1) = 0;
end