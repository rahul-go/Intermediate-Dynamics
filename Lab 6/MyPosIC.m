function [E] = MyPosIC(r, t_2, x)

%% Solved Values
% See attached file for hand calculations.

% Easy access to...
t_3 = x(1);                     % Angular position of link 3
t_4 = x(2);                     % Angular position of link 4

% Find the Error
e_x = r(1) + r(4)*cos(t_4) - r(2)*cos(t_2) - r(3)*cos(t_3);
e_y = r(4)*sin(t_4) - r(2)*sin(t_2) - r(3)*sin(t_3);
E = hypot(e_x, e_y);

end