%% BowlingBallEOM Usage and Description
% TODO



function [xdot] = BowlingBallEOM(x)

    % Given Values
    w = 15;                     % Weight of bowling ball (lb)
    r = (8.5/2)/12;             % Radius of bowling ball (ft)
    mu_k = 0.12;                % Kinetic friction coefficient (unitless)
    mu_s = 0.14;                % Static friction coefficient (unitless)
    g = 32.174;                 % Gravitational constant (in/s^2)
    
    % Solved Values
    m = w/g;                    % Mass of bowling ball (slug)
    
    
    
    v_x = x(1);                 % Easy access to v_x
    v_y = x(2);                 % Easy access to v_y
    omega_x = x(3);             % Easy access to omega_x
    omega_y = x(4);             % Easy access to omega_y
    
    N = m*g;                    % Normal force (lb)
    
    vc_x = v_x - omega_y*r;     % Velocity[x] (m/s)
    vc_y = v_y + omega_x*r;     % Velocity[y] (m/s)
    theta = atan2(vc_y,vc_x);   % Theta (deg)

    F_x = -mu_k*N*cosd(theta);  % Friction force[x] (lb)
    F_y = -mu_k*N*sind(theta);  % Friction force[y] (lb)

    I = 2/5*m*r^2;              % Moment of inertia (slug*ft^2)

    xdot = [F_x/m;
            F_y/m;
            F_y*r/I;
            F_x*r/I;
            v_x;
            v_y];

end