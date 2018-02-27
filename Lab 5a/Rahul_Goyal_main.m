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

% Set Initial theta_2
% t2_0 = deg2rad(140);            % Angular position of TODO (rad)
t2_0 = deg2rad(0);              % Angular position of TODO (rad)

% Position Initial Conditions
r3_0 = hypot(r_1+r_2*cos(t2_0), r_2*sin(t2_0));
                                % TODO (m) [Pythagorean Theorem]
t3_0 = asin(sin(pi-t2_0)/r3_0 * r_2);
                                % Theta_3 Initial (rad) [Law of Sines]
x2_0 = r_1+r_2*cos(t2_0)/2;     % TODO (m)
y2_0 = r_2*sin(t2_0)/2;         % TODO (m)
x3_0 = r_1+r_2*cos(t2_0) - l_ab*cos(t3_0)/2;
                                % TODO (m)
y3_0 = r_2*sin(t2_0) - l_ab*sin(t3_0)/2;
                                % TODO (m)

% Position Initial Conditions Matrix
% x_0 = [r3_0, t3_0];
x_0 = [r3_0, t3_0, x2_0, y2_0, x3_0, y3_0];

% Velocity Initial Conditions
rdot3_0 = 0;                    % TODO (m/s) FIX
tdot3_0 = r_2*tdot_2 / r3_0;    % TODO (rad/s) FIX
xdot2_0 = tdot_2*r_2*sin(t2_0) / 2;
                                % TODO (m/s)
ydot2_0 = tdot_2*r_2*cos(t2_0) / 2;
                                % TODO (m/s)
xdot3_0 = tdot3_0*r3_0*sin(t3_0) / 2;
                                % TODO (m/s)
ydot3_0 = tdot3_0*r3_0*cos(t3_0) / 2;
                                % TODO (m/s)

% Velocity Initial Conditions Matrix
% v_0 = [rdot3_0, tdot3_0];
v_0 = [rdot3_0, tdot3_0, xdot2_0, ydot2_0, xdot3_0, ydot3_0];



%% Simulate the Slider-Crank Using Simulink
% TODO
sim('untitled.slx');



%% Plotting Data Setup
% TODO

% TODO
r1_x = [0, r_1];
r1_y = [0, 0];



%% Tip-to-Tail Animation
% TODO

% for t = 1:length(tout)
% 
%     % TODO
%     t_2 = t2_0 + tdot_2 * tout(t);
% 
%     % TODO
%     r2_x = [r1_x(end), r1_x(end) + r_2*cos(t_2)];
%     r2_y = [r1_y(end), r1_y(end) + r_2*sin(t_2)];
%     r3_x = [r2_x(end), r1_x(1)];
%     r3_y = [r2_y(end), r2_y(1)];
% 
%     % Plot the links
%     plot(r1_x, r1_y, r2_x, r2_y, r3_x, r3_y, 'LineWidth', 2);
%     % Keeps the frame consistent
%     axis equal;
%     axis([-0.2, 0.4, -0.1, 0.1]);
%     
%     % Calculate the time step and pause accordingly
%     if t ~= length(tout)            % Prevent index error
%         % Calculate the time step (s)
%         t_step = tout(t+1) - tout(t);
%         pause(t_step);              % Assume negligible processing time
%     end
% 
% end



%% Plot of R3 vs. theta_2
% TODO
plot(tdot_2*tout, vout(:, 1));



%% Simulation Animation
% TODO

for t = 1:length(tout)

    % TODO
    t_2 = t2_0 + tdot_2 * tout(t);
    
    % TODO
    r_3 = xout(t, 1);               % Easy access to TODO
    t_3 = xout(t, 2);               % Easy access to TODO
    x_2 = xout(:, 3);               % Easy access to TODO
    y_2 = xout(:, 4);               % Easy access to TODO
    x_3 = xout(:, 5);               % Easy access to TODO
    y_3 = xout(:, 6);               % Easy access to TODO

    % TODO
    r2_x = [r1_x(end), r1_x(end) + r_2*cos(t_2)];
    r2_y = [r1_y(end), r1_y(end) + r_2*sin(t_2)];
    lab_x = [r2_x(end), r2_x(end) - l_ab*cos(t_3)];
    lab_y = [r2_y(end), r2_y(end) - l_ab*sin(t_3)];

    % Plot
    plot(r1_x, r1_y, ...            % Link 1
         r2_x, r2_y, ...            % Link 2
         lab_x, lab_y, ...          % Link AB
         x_2(1:t), y_2(1:t), ...    % TODO
         x_3(1:t), y_3(1:t), ...    % TODO
         'LineWidth', 2);           % Line Properties

    % Keeps the frame consistent
    axis equal;
    axis([-0.2, 0.4, -0.1, 0.1]);
    
    % Calculate the time step and pause accordingly
    if t ~= length(tout)            % Prevent index error
        % Calculate the time step (s)
        t_step = tout(t+1) - tout(t);
        pause(t_step);              % Assume negligible processing time
    end

end