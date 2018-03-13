function [E] = MyPosIC(t_2, x)

%% Set Values
% The following is used to easily change the lengths of the vectors. (A
% Grashof mechanism has the constraint R1 + R2 <= R3 + R4).
r_1 = 2;                        % Length of link 1 (m)
r_2 = 4;                        % Length of link 2 (m)
r_3 = 6;                        % Length of link 3 (m)
r_4 = 8;                        % Length of link 4 (m)



%% Solved Values
% See attached file for hand calculations.

% Easy access to...
t_3 = x(1);                     % Angular position of link 3
t_4 = x(2);                     % Angular position of link 4

% Find the Error
e_x = r_1 + r_4*cos(t_4) - r_2*cos(t_2) - r_3*cos(t_3);
e_y = r_4*sin(t_4) - r_2*sin(t_2) - r_3*sin(t_3);
E = hypot(e_x, e_y);



end