%% link_solver Usage and Description
% TODO

function [x] = link_solver(u)

%% Given Values
% The following assigns values given by the problem statement to variables.

% Given Values
r_1 = 240/1000;                 % TODO (m)
r_2 = 80/1000;                  % TODO (m)
l_ab = 350/1000;                % TODO (m)



%% Easy Access
% TODO
% t_2 = u(1);
% tdot_2 = u(2);
% r_3 = u(3);
% t_3 = u(4);
% rdot_3 = u(5);
% tdot_3 = u(6);

t_2 = u(1);
tdot_2 = u(2);
r_3 = u(3);
t_3 = u(4);
x_2 = u(5);
y_2 = u(6);
x_3 = u(7);
y_3 = u(8);
rdot_3 = u(9);
tdot_3 = u(10);
xdot_2 = u(11);
ydot_2 = u(12);
xdot_3 = u(13);
ydot_3 = u(14);



%% Solved Values
% The following assigns values derived and/or solved from the given values
% to variables. See the attached file for hand calculations.

% A = [cos(t_3), -r_3*sin(t_3);
%      sin(t_3), r_3*cos(t_3)];
 
A = [cos(t_3), -r_3*sin(t_3), 0, 0, 0, 0;
     sin(t_3), r_3*cos(t_3), 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     -cos(t_3), (r_3-l_ab/2)*sin(t_3), 0, 0, 1, 0;
     -sin(t_3), -(r_3-l_ab/2)*cos(t_3), 0, 0, 0, 1];

% B = [2*rdot_3*tdot_3*sin(t_3) + r_3*tdot_3^2*cos(t_3) - r_2*tdot_2^2*cos(t_2);
%      -2*rdot_3*tdot_3*cos(t_3) + r_3*tdot_3^2*sin(t_3) - r_2*tdot_2^2*sin(t_2)];

B = [2*rdot_3*tdot_3*sin(t_3) + r_3*tdot_3^2*cos(t_3) - r_2*tdot_2^2*cos(t_2);
     -2*rdot_3*tdot_3*cos(t_3) + r_3*tdot_3^2*sin(t_3) - r_2*tdot_2^2*sin(t_2);
     -1/2*r_2*tdot_2^2*cos(t_2);
     -1/2*r_2*tdot_2^2*sin(t_2);
     -2*rdot_3*tdot_3*sin(t_3) - r_3*tdot_3^2*cos(t_3) + l_ab/2*tdot_3^2*cos(t_3);
     2*rdot_3*tdot_3*cos(t_3) - r_3*tdot_3^2*sin(t_3) + l_ab/2*tdot_3^2*sin(t_3)];



%% Solve for x
% TODO
x = A \ B;

end