%% Rahul_Goyal_main Usage and Description
% ME 326 Winter 2018 - Laboratory Assignment #5
%
% *Author:* RAHUL GOYAL
%
% California Polytechnic State University, San Luis Obispo, CA
%
% *Date Created:* February 13, 2018
%
% *Date Modified:* February 27, 2018
%
% *Description:*
% TODO
% 
% *Required Files:*
%
% * TODO
%
% *Still To Do:*
%
% * Start!



%% Problem Statement
% TODO



%% Reset
% The following was used while debugging.

close all;
clear all;
clc;



%% Given Values
% The following assigns values given by the problem statement to variables.

% Given Values
tdot_2 = -3;                    % Angular velocity of TODO (rad/s)
l_ab = 350/1000;                % TODO (m)
r_1 = 240/1000;                 % TODO (m)
r_2 = 80/1000;                  % TODO (m)
theta2_stop = -4*pi;            % TODO (rad)



%% Initial Conditions
% The following sets the initial conditions of the slider-crank. See the
% attached file for hand calculations.

% TODO
% t2_0 = deg2rad(140);            % Angular position of TODO (rad)
t2_0 = deg2rad(0);              % Angular position of TODO (rad)

% Position Initial Conditions
r3_0 = hypot(r_1+r_2*cos(t2_0), r_2*sin(t2_0));
                                % TODO (m) [Pythagorean Theorem]
t3_0 = asin(sin(pi-t2_0)/r3_0 * r_2);
                                % Theta_3 Initial (rad) [Law of Sines]
x_0 = [r3_0, t3_0];
% x_0 = [r3_0, t3_0, x2_0, y2_0, x3_0, y3_0];

% Velocity Initial Conditions
rdot3_0 = 0;                    % TODO (m/s) FIX
tdot3_0 = r_2*tdot_2 / r3_0;    % TODO (rad/s) FIX
v_0 = [rdot3_0, tdot3_0];
% v_0 = [rdot3_0, tdot3_0, xdot2_0, ydot2_0, xdot3_0, ydot3_0];



%% Simulate the Slider-Crank Using Simulink
% TODO
sim('untitled.slx');



%% Tip-to-Tail Animation
% TODO

% TODO
r1_x = [0, r_1];
r1_y = [0, 0];

for t = 1:length(tout)

    % TODO
    t_2 = t2_0 + tdot_2 * tout(t);
    
%     % TODO
%     r_3 = x_out(t, 1);          % Easy access to TODO
%     t_3 = x_out(t, 2);          % Easy access to TODO

    % TODO
    r2_x = [r1_x(end), r1_x(end) + r_2*cos(t_2)];
    r2_y = [r1_y(end), r1_y(end) + r_2*sin(t_2)];
    r3_x = [r2_x(end), r1_x(1)];
    r3_y = [r2_y(end), r2_y(1)];

    % Plot the links
    plot(r1_x, r1_y, r2_x, r2_y, r3_x, r3_y, 'LineWidth', 2);
    % Keeps the frame consistent
    axis equal;
    axis([-0.2, 0.4, -0.1, 0.1]);
%     pause(0.1);                 % TODO TEMP

end



%% Plot of R3 vs. theta_2
% TODO
plot(tdot_2*tout, v_out(:, 1));
pause(5);



%% Simulation Animation
% TODO

% TODO
r1_x = [0, r_1];
r1_y = [0, 0];

for t = 1:length(tout)

    % TODO
    t_2 = t2_0 + tdot_2 * tout(t);
    
    % TODO
    r_3 = x_out(t, 1);          % Easy access to TODO
    t_3 = x_out(t, 2);          % Easy access to TODO

    % TODO
    r2_x = [r1_x(end), r1_x(end) + r_2*cos(t_2)];
    r2_y = [r1_y(end), r1_y(end) + r_2*sin(t_2)];
    r3_x = [r2_x(end), r2_x(end) - l_ab*cos(t_3)];
    r3_y = [r2_y(end), r2_y(end) - l_ab*sin(t_3)];

    % Plot the links
    plot(r1_x, r1_y, r2_x, r2_y, r3_x, r3_y, 'LineWidth', 2);
    % Keeps the frame consistent
    axis equal;
    axis([-0.2, 0.4, -0.1, 0.1]);
    pause(0.2);                 % TODO TEMP

end