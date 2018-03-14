function [E] = MyPosIC(r, t_2, x)

%% Declare Global Variables
% The following declares global variables.
global image t_step;



%% Solved Values
% See attached file for hand calculations.

% Easy access to...
t_3 = x(1);                     % Angular position of link 3
t_4 = x(2);                     % Angular position of link 4

% Find the Error
e_x = r(1) + r(4)*cos(t_4) - r(2)*cos(t_2) - r(3)*cos(t_3);
e_y = r(4)*sin(t_4) - r(2)*sin(t_2) - r(3)*sin(t_3);
E = hypot(e_x, e_y);



%% Error Animation
% The following animates the error in the angular positions of the links of
% the crank-rocker.

% Cartesian Coordinates of Link 1
r1_x = [0, r(1)];
r1_y = [0, 0];
% Cartesian Coordinates of Link 2
r2_x = [0, r(2)*cos(t_2)];
r2_y = [0, r(2)*sin(t_2)];
% Cartesian Coordinates of Link 3
r3_x = [r2_x(end), r2_x(end) + r(3)*cos(t_3)];
r3_y = [r2_y(end), r2_y(end) + r(3)*sin(t_3)];
% Cartesian Coordinates of Link 4
r4_x = [r1_x(end), r1_x(end) + r(4)*cos(t_4)];
r4_y = [r1_y(end), r1_y(end) + r(4)*sin(t_4)];

% Plot the links
plot(r1_x, r1_y, ...            % Link 1
     r2_x, r2_y, ...            % Link 2
     r3_x, r3_y, ...            % Link 3
     r4_x, r4_y, ...            % Link 4
     'LineWidth', 2);           % Line Properties
% Keep the frame consistent
axis equal;
axis([-r(2)-0.5, r(1)+r(4)+5, -r(4)-0.5, r(4)+0.5]);
% Store the time step for later use
t_step(length(t_step)+1) = 0.0001;
% pause(t_step(end));             % Pause for humans

% Convert the plot frame to an image and store for later use
image{length(image)+1} = frame2im(getframe(1));

end