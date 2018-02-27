%% link_solver Usage and Description
% TODO

function [x] = link_solver(u)

%% Given Values
% The following assigns values given by the problem statement to variables.

% Given Values
l_ab = 350/1000;                % Length of link AB (m)
r_1 = 240/1000;                 % Length of vector R1 (m)
r_2 = 80/1000;                  % Length of vector R2 (m)



%% Solved Values
% The following assigns values derived and/or solved from the given values
% to variables. See the attached file for hand calculations.

% Easy Access to...
t_2 = u(1);                     % Angular position of link OA (rad)
tdot_2 = u(2);                  % Angular velocity of link OA (rad/s)
r_3 = u(3);                     % Length of vector R3 (m)
t_3 = u(4);                     % Angular position of link AB (rad)
rdot_3 = u(9);                  % Velocity of vector R3 (m/s)
tdot_3 = u(10);                 % Angular velocity of link AB (rad/s)

A = [cos(t_3), -r_3*sin(t_3), 0, 0, 0, 0;
     sin(t_3), r_3*cos(t_3), 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     -cos(t_3), (r_3-l_ab/2)*sin(t_3), 0, 0, 1, 0;
     -sin(t_3), -(r_3-l_ab/2)*cos(t_3), 0, 0, 0, 1];

b = [2*rdot_3*tdot_3*sin(t_3) + r_3*tdot_3^2*cos(t_3) - r_2*tdot_2^2*cos(t_2);
     -2*rdot_3*tdot_3*cos(t_3) + r_3*tdot_3^2*sin(t_3) - r_2*tdot_2^2*sin(t_2);
     -1/2*r_2*tdot_2^2*cos(t_2);
     -1/2*r_2*tdot_2^2*sin(t_2);
     -2*rdot_3*tdot_3*sin(t_3) - r_3*tdot_3^2*cos(t_3) + l_ab/2*tdot_3^2*cos(t_3);
     2*rdot_3*tdot_3*cos(t_3) - r_3*tdot_3^2*sin(t_3) + l_ab/2*tdot_3^2*sin(t_3)];



%% Solve for x
% TODO
x = A \ b;

end