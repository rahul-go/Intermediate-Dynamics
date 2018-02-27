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
% * Simulator.slx - TODO
% * link_solver.m - TODO
%
% *Still To Do:*
%
% * Start!



%% Problem Statement
% Crank OA rotates clockwise at the constant rate theta-dot = 3 rad/s. The
% connecting link AB passes through the pivoted collar at C. Note: lAB =
% 350 mm. Develop a complete kinematic model of the slider-crank mechanism.



%% Reset
% The following was used while debugging.

close all;
clear all;
clc;



%% Given Values
% The following assigns values given by the problem statement to variables.

% Given Values
tdot_2 = -3;                    % Angular velocity of link OA (rad/s)
l_ab = 350/1000;                % Length of link AB (m)
r_1 = 240/1000;                 % Length of vector R1 (m)
r_2 = 80/1000;                  % Length of vector R2 (m)
theta2_stop = -4*pi;            % Theta_2 final (rad)

% Set Values
t2_0 = deg2rad(0);              % Angular position initial of link OA (rad)
% t2_0 = deg2rad(140);            % Angular position initial of link OA (rad)



%% Initial Conditions
% The following sets the initial conditions of the slider-crank. See the
% attached file for hand calculations.

% Position Initial Conditions
% Length of vector R2 initial (m) [Pythagorean Theorem]
r3_0 = hypot(r_1+r_2*cos(t2_0), r_2*sin(t2_0));
% Angular position of link OA initial (rad) [Law of Sines]
t3_0 = asin(sin(pi-t2_0)/r3_0 * r_2);
% COM[x] initial of link OA (m)
x2_0 = r_1+r_2*cos(t2_0)/2;
% COM[y] initial of link OA (m)
y2_0 = r_2*sin(t2_0)/2;
% COM[x] initial of link AB (m)
x3_0 = r_1+r_2*cos(t2_0) - l_ab*cos(t3_0)/2;
% COM[y] initial of link AB (m)
y3_0 = r_2*sin(t2_0) - l_ab*sin(t3_0)/2;

% Position Initial Conditions Matrix
x_0 = [r3_0, t3_0, x2_0, y2_0, x3_0, y3_0];



% Velocity Initial Conditions
% Velocity of vector R2 initial (m/s) FIX
rdot3_0 = 0;
% Angular velocity of link OA initial (rad/s) FIX
tdot3_0 = r_2*tdot_2 / r3_0;
% Velocity_G[x] initial of link OA (m/s)
xdot2_0 = tdot_2*r_2*sin(t2_0) / 2;
% Velocity_G[y] initial of link OA (m/s)
ydot2_0 = tdot_2*r_2*cos(t2_0) / 2;
% Velocity_G[x] initial of link AB (m/s)
xdot3_0 = tdot3_0*r3_0*sin(t3_0) / 2;
% Velocity_G[x] initial of link AB (m/s)
ydot3_0 = tdot3_0*r3_0*cos(t3_0) / 2;

% Velocity Initial Conditions Matrix
v_0 = [rdot3_0, tdot3_0, xdot2_0, ydot2_0, xdot3_0, ydot3_0];



%% Simulate the Slider-Crank Using Simulink
% TODO
sim('Simulator.slx');



%% Plotting Data Setup
% TODO
r1_x = [0, r_1];
r1_y = [0, 0];



%% Tip-to-Tail Animation
% TODO

for t = 1:length(tout)

    t_2 = t2_0 + tdot_2 * tout(t);  % Angular position of link OA (m)

    % Cartesian Coordinates of Link OA, AB (tip-to-tail)
    r2_x = [r1_x(end), r1_x(end) + r_2*cos(t_2)];
    r2_y = [r1_y(end), r1_y(end) + r_2*sin(t_2)];
    r3_x = [r2_x(end), r1_x(1)];
    r3_y = [r2_y(end), r2_y(1)];

    % Plot the vector links
    plot(r1_x, r1_y, r2_x, r2_y, r3_x, r3_y, 'LineWidth', 2);
    % Keep the frame consistent
    axis equal;
    axis([-0.2, 0.4, -0.1, 0.1]);
    
    % Calculate the time step and pause accordingly
    if t ~= length(tout)            % Prevent index error
        % Calculate the time step (s)
        t_step = tout(t+1) - tout(t);
        pause(t_step);              % Assume negligible processing time
    end

end



%% Plot of R3 vs. Theta_2
% TODO
plot(tdot_2*tout, vout(:, 1));



%% Simulation Animation
% TODO

for t = 1:length(tout)

    t_2 = t2_0 + tdot_2*tout(t);    % Angular position of link OA (m)
    
    % Easy Access to...
    r_3 = xout(t, 1);               % Length of vector R3 (m)
    t_3 = xout(t, 2);               % Angular position of link AB (rad)
    x_2 = xout(:, 3);               % COM[x] of link OA (m)
    y_2 = xout(:, 4);               % COM[y] of link OA (m)
    x_3 = xout(:, 5);               % COM[x] of link AB (m)
    y_3 = xout(:, 6);               % COM[y] of link AB (m)

    % Cartesian Coordinates of Link OA, AB
    r2_x = [r1_x(end), r1_x(end) + r_2*cos(t_2)];
    r2_y = [r1_y(end), r1_y(end) + r_2*sin(t_2)];
    lab_x = [r2_x(end), r2_x(end) - l_ab*cos(t_3)];
    lab_y = [r2_y(end), r2_y(end) - l_ab*sin(t_3)];

    % Plot the links, COMs, COM paths
    plot(r1_x, r1_y, ...            % Link 1
         r2_x, r2_y, ...            % Link 2
         lab_x, lab_y, ...          % Link AB
         x_2(1:t), y_2(1:t), ...    % Path of link OA COM
         x_3(1:t), y_3(1:t), ...    % Path of link AB COM
         'LineWidth', 2);           % Line Properties

    % Keep the frame consistent
    axis equal;
    axis([-0.2, 0.4, -0.1, 0.1]);
    
    % Calculate the time step and pause accordingly
    if t ~= length(tout)            % Prevent index error
        % Calculate the time step (s)
        t_step = tout(t+1) - tout(t);
        pause(t_step);              % Assume negligible processing time
    end

end