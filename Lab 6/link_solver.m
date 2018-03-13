%% link_solver Usage and Description
% This function represents the accelerations of the simulation. It returns
% x with an input of u.

function [x] = link_solver(u)

%% Set Values
% The following is used to easily change the lengths and masses of the
% links. (A Grashof mechanism has the constraint R1 + R2 <= R3 + R4).
r = [2, 3, 4, 5];               % Length of links 1, 2, 3, 4 (m)
m_2 = 1;                        % Mass of link 2 (kg)
m_3 = 1;                        % Mass of link 3 (kg)
m_4 = 1;                        % Mass of link 4 (kg)



%% Easy Access
% The following sets up easy access to variables.
t_2 = u(1);                     % TODO
tdot_2 = u(2);                  % TODO
t_3 = u(3);                     % TODO
t_4 = u(4);                     % TODO
tdot_3 = u(11);                 % TODO
tdot_4 = u(12);                 % TODO

% Simplicity and Compactness of Notation
c_2 = cos(t_2);
s_2 = sin(t_2);
c_3 = cos(t_3);
s_3 = sin(t_3);
c_4 = cos(t_4);
s_4 = sin(t_4);



%% Solved Values
% The following assigns values derived and/or solved from the given values
% to variables. See the attached file for hand calculations.

I_3 = 1/12*m_3*r(3)^2;          % Moment of inertia of link 3 (kg*m^2)
I_4 = 1/12*m_4*r(4)^2;          % Moment of inertia of link 4 (kg*m^2)

A = [-r(3)*s_3, r(4)*s_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     r(3)*c_3, -r(4)*c_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     r(3)/2*s_3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     -r(3)/2*c_3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, r(4)/2*s_4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, -r(4)/2*c_4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, m_2, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0;
     0, 0, 0, m_2, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, r(2)/2*s_2, -r(2)/2*c_2, -r(2)/2*s_2, r(2)/2*c_2, 0, 0, 0, 0, -1;
     0, 0, 0, 0, m_3, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0;
     0, 0, 0, 0, 0, m_3, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0;
     I_3, 0, 0, 0, 0, 0, 0, 0, r(3)/2*s_3, -r(3)/2*c_3, 0, 0, 0, 0, r(3)/2*s_3, -r(3)/2*c_3, 0;
     0, 0, 0, 0, 0, 0, m_4, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0;
     0, 0, 0, 0, 0, 0, 0, m_4, 0, 0, 0, 0, 0, -1, 0, 1, 0;
     0, I_4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -r(4)/2*s_4, r(4)/2*c_3, -r(4)/2*s_4, r(4)/2*c_4, 0];

b = [-r(4)*c_4*tdot_4^2 + r(2)*c_2*tdot_2^2 + r(3)*c_3*tdot_3^2;
     -r(4)*s_4*tdot_4^2 + r(2)*s_2*tdot_2^2 + r(3)*s_3*tdot_3^2;
     -r(2)/2*c_2*tdot_2^2;
     -r(2)/2*s_2*tdot_2^2;
     -r(2)*c_2*tdot_2^2 - r(3)/2*c_3*tdot_3^2;
     -r(2)*s_2*tdot_2^2 - r(3)/2*s_3*tdot_3^2;
     -r(4)/2*c_4*tdot_4^2;
     -r(4)/2*s_4*tdot_4^2;
     0;
     0;
     0;
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