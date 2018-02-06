%% CarEOM Usage and Description
% This function represents the equations of motion for the simulation. It
% returns xdot with inputs of time, x, the A and B state-space matrices,
% the slope data, and the horizontal velocity of the car. See the attached
% sheet for hand calculations. TODO UPDATE
%
% To find u (in this case, the vertical component of the velocity of
% the car) corresponding to the time and x value, this function solves for
% the corresponding distance value, from which it interpolates the
% corresponding slope value, from which it solves for the corresponding u
% value, the vertical component of the velocity of the car. Then, knowing
% the x and u values, the script solves y = A*x + B*u using matrix
% multiplication.
% 
% To unstick the tires, this function TODO

function [xdot] = CarEOM(t, x, A, B, slope, v_x, w_c, w_u, k_s, b_s, k_t, m_u)
    
    x_r = v_x*t;                % Distance of car (ft)
    % Interpolate slope (m) from distance (in/ft)
    m_r = interp1(slope(:, 1), slope(:, 2), x_r);
    u = v_x*m_r;                % Vertical velocity (yr_dot) of car (ft/s)
    
    
    
    xdot = A*x + B*u;
    
    
    
    % Unstuck tires
    y_r = x(1);                 % Easy access to y_r
    y_u = x(2);                 % Easy access to y_u
    
    f_kt = k_t*(y_u-y_r);       % Force in spring k_t (stuck)
    
    if f_kt > (w_c+w_u)
        
        f_kt = w_c+w_u;             % Force in spring k_t (unstuck)
        
        y_c = x(3);                 % Easy access to y_c
        yu_dot = x(4);              % Easy access to yu_dot
        yc_dot = x(5);              % Easy access to yc_dot
        
        f_ks = k_s*(y_c-y_u);       % Force in spring k_s
        f_bs = b_s*(yc_dot-yu_dot); % Force in damper b_s
        
        yu_2dot = (f_ks + f_bs - f_kt)/m_u;
        
        xdot(4) = yu_2dot;
        
    end
    
end