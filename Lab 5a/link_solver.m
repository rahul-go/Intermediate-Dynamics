%% link_solver Usage and Description
% TODO

function [x] = link_solver(u)

%% Given Values
% The following assigns values given by the problem statement to variables.

% Given Values
r_1 = 240/1000;                 % TODO (m)
r_2 = 80/1000;                  % TODO (m)



%% Easy Access
% TODO
t_2 = u(1);
tdot_2 = u(2);
r_3 = u(3);
t_3 = u(4);
rdot_3 = u(5);
tdot_3 = u(6);



%% Solved Values
% The following assigns values derived and/or solved from the given values
% to variables. See the attached sheet for hand calculations.

A = [cos(t_3), -r_3*sin(t_3);
     sin(t_3), r_3*cos(t_3)];

B = [2*rdot_3*tdot_3*sin(t_3) + r_3*tdot_3^2*cos(t_3) - r_2*tdot_2^2*cos(t_2);
     -2*rdot_3*tdot_3*cos(t_3) + r_3*tdot_3^2*sin(t_3) - r_2*tdot_2^2*sin(t_3)];



%% Solve for x
% TODO
x = A \ B;

end