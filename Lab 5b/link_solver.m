%% link_solver Usage and Description
% This function represents the accelerations of the simulation. It returns
% x with an input of u.

function [x] = link_solver(u)

%% Given Values
% The following assigns values given by the problem statement to variables.
l_ab = 350/1000;                % Length of link AB (m)
r_1 = 240/1000;                 % Length of vector R1 (m)
r_2 = 80/1000;                  % Length of vector R2 (m)
m_2 = 10;                       % Mass of link OA (kg)
m_3 = 15;                       % Mass of link AB (kg)



%% Solved Values
% The following assigns values derived and/or solved from the given values
% to variables. See the attached file for hand calculations.

% Easy access to...
t_2 = u(1);                     % Angular position of link OA (rad)
tdot_2 = u(2);                  % Angular velocity of link OA (rad/s)
r_3 = u(3);                     % Length of vector R3 (m)
t_3 = u(4);                     % Angular position of link AB (rad)
rdot_3 = u(9);                  % Velocity of vector R3 (m/s)
tdot_3 = u(10);                 % Angular velocity of link AB (rad/s)

I_3 = 1/12*m_3*l_ab^2;          % Moment of inertia of link AB (kg*m^2)

A = [cos(t_3), -r_3*sin(t_3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     sin(t_3), r_3*cos(t_3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
     -cos(t_3), (r_3-l_ab/2)*sin(t_3), 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
     -sin(t_3), -(r_3-l_ab/2)*cos(t_3), 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
     0, 0, m_2, 0, 0, 0, 1, 0, -1, 0, 0, 0;
     0, 0, 0, m_2, 0, 0, 0, 1, 0, -1, 0, 0;
     0, 0, 0, 0, 0, 0, -r_2*sin(t_2)/2, r_2*cos(t_2)/2, -r_2*sin(t_2)/2, r_2*cos(t_2)/2, 0, 1;
     0, 0, 0, 0, m_3, 0, 0, 0, 1, 0, sin(t_3), 0;
     0, 0, 0, 0, 0, m_3, 0, 0, 0, 1, -cos(t_3), 0;
     0, I_3, 0, 0, 0, 0, 0, 0, -l_ab*sin(t_3)/2, l_ab*cos(t_3)/2, r_3-l_ab/2, 0];

b = [2*rdot_3*tdot_3*sin(t_3) + r_3*tdot_3^2*cos(t_3) - r_2*tdot_2^2*cos(t_2);
     -2*rdot_3*tdot_3*cos(t_3) + r_3*tdot_3^2*sin(t_3) - r_2*tdot_2^2*sin(t_2);
     -1/2*r_2*tdot_2^2*cos(t_2);
     -1/2*r_2*tdot_2^2*sin(t_2);
     -2*rdot_3*tdot_3*sin(t_3) - r_3*tdot_3^2*cos(t_3) + l_ab/2*tdot_3^2*cos(t_3);
     2*rdot_3*tdot_3*cos(t_3) - r_3*tdot_3^2*sin(t_3) + l_ab/2*tdot_3^2*sin(t_3);
     0;
     0;
     0;
     0;
     0;
     0];



%% Solve for x
% Solve for x using mldivide.
x = A \ b;

end