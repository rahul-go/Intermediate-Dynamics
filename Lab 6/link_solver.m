%% link_solver Usage and Description
% This function represents the accelerations of the simulation. It returns
% x with an input of u.

function [x] = link_solver(u)

%% Given Values
% The following assigns values given by the problem statement to variables.

% Given Values
r_1 = 2;                        % TODO
r_2 = 3;                        % TODO
r_3 = 4;                        % TODO
r_4 = 5;                        % TODO



%% Easy Access
% The following sets up easy access to variables.
t_2 = u(1);                     % TODO
tdot_2 = u(2);                  % TODO
t_3 = u(3);                     % TODO
t_4 = u(4);                     % TODO
tdot_3 = u(5);                  % TODO
tdot_4 = u(6);                  % TODO



%% Solved Values
% The following assigns values derived and/or solved from the given values
% to variables. See the attached file for hand calculations.

c_2 = cos(t_2);
s_2 = sin(t_2);
c_3 = cos(t_3);
s_3 = sin(t_3);
c_4 = cos(t_4);
s_4 = sin(t_4);

A = [-r_3*s_3, r_4*s_4;
     r_3*c_3, -r_4*c_4];

b = [-r_4*c_4*tdot_4^2 + r_2*c_2*tdot_2^2 + r_3*c_3*tdot_3^2;
     -r_4*s_4*tdot_4^2 + r_2*s_2*tdot_2^2 + r_3*s_3*tdot_3^2];



%% Solve for x
% Solve for x using mldivide.
x = A \ b;

end