%% Rahul_Goyal_main Usage and Description
% ME 326 Winter 2018 - Laboratory Assignment #6
%
% *Author:* RAHUL GOYAL
%
% California Polytechnic State University, San Luis Obispo, CA
%
% *Date Created:* March 06, 2018
%
% *Date Modified:* March 13, 2018
%
% *Description:*
% This script simulates the motion of a crank-rocker. TODO
% 
% *Required Files:*
%
% * MyPosIC.m - TODO
% * Simulator.slx - This file uses Simulink to double integrate part of a
% MATLAB Function Block which describes the accelerations and forces of the
% simulation (only the accelerations are integrated). It outputs the
% positions as xout, the velocities as vout, the accelerations as aout, the
% forces as Fout, and the times as tout with inputs of the MATLAB function,
% initial conditions, and final condition.
% * link_solver.m - This file contains a function that represents the
% accelerations and forces of the simulation. It returns x with an input of
% u.
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



%% Set Values
% The following is used to easily change the lengths of the vectors and the
% angular velocity of link 2. (A Grashof mechanism has the constraint R1 +
% R2 <= R3 + R4).
r = [2, 3, 4, 5];               % Length of links 1, 2, 3, 4 (m)
tdot_2 = 1;                     % Angular velocity of link 2 (rad/s)



%% Given Values
% The following assigns values given by the problem statement to variables.
t2_f = -4*pi;                   % TODO (rad)



%% Initial Conditions
% The following sets the initial conditions of the crank-rocker.

% Position Initial Conditions
t2_0 = deg2rad(30);             % Angular position initial of link 2 (rad)
t3_0 = 0;                       % Angular position initial of link 3 (rad) [GUESS]
t4_0 = 0;                       % Angular position initial of link 4 (rad) [GUESS]
% TODO
minimize = @(x) MyPosIC(r, t2_0, x);

% Position Initial Conditions Matrix
x_0 = fminsearch(minimize, [t3_0, t4_0]);

% Easy access to...
t3_0 = x_0(1);                  % Angular position initial of link 3 (rad)
t4_0 = x_0(2);                  % Angular position initial of link 4 (rad)

% Simplicity and Compactness of Notation
c2_0 = cos(t2_0);
s2_0 = sin(t2_0);
c3_0 = cos(t3_0);
s3_0 = sin(t3_0);
c4_0 = cos(t4_0);
s4_0 = sin(t4_0);

% Velocity Initial Conditions
A = [-r(3)*s3_0, r(4)*s4_0
     r(3)*c3_0, -r(4)*c4_0];
b = [r(2)*s2_0*tdot_2
     -r(2)*c2_0*tdot_2];

% Velocity Initial Conditions Matrix
v_0 = A \ b;



%% Simulate the Crank-Rocker Using Simulink
% The following calls the Simulink file Simulator.slx, which outputs the
% positions as xout, the velocities as vout, the accelerations as aout, the
% forces as Fout, and the times as tout with with link_solver.m as the
% input for the MATLAB Fuction, tdot_2, t2_0, v_0, and x_0 as the inputs
% for the initial conditions, and t2_f as the input for the final
% conditions.
sim('Simulator.slx');



%% Simulation Animation
% The following animates the crank-rocker by using the simulation data.

% Cartesian Coordinates of Link 1
r1_x = [0, r(1)];
r1_y = [0, 0];

% Easy access to...
t_2 = t2_0 + tdot_2*tout;       % Angular position of link 2 (rad)
t_3 = xout(:, 1);               % Angular position of link 3 (rad)
t_4 = xout(:, 2);               % Angular position of link 4 (rad)

for t = 1:length(tout)
    
    % Cartesian Coordinates of Link 2
    r2_x = [0, r(2)*cos(t_2(t))];
    r2_y = [0, r(2)*sin(t_2(t))];
    % Cartesian Coordinates of Link 3
    r3_x = [r2_x(end), r2_x(end) + r(3)*cos(t_3(t))];
    r3_y = [r2_y(end), r2_y(end) + r(3)*sin(t_3(t))];
    % Cartesian Coordinates of Link 4
    r4_x = [r1_x(end), r1_x(end) + r(4)*cos(t_4(t))];
    r4_y = [r1_y(end), r1_y(end) + r(4)*sin(t_4(t))];
    
    % TODO
    plot(r1_x, r1_y, ...            % Link 1
         r2_x, r2_y, ...            % Link 2
         r3_x, r3_y, ...            % Link 3
         r4_x, r4_y, ...            % Link 4
         'LineWidth', 2);           % Line Properties
    % Keep the frame consistent
    axis equal;
    axis([-r(2)-0.5, r(1)+r(4)+0.5, -r(4)-0.5, r(4)+0.5]);
    
	% Calculate the time step and pause accordingly
    if t ~= length(tout)            % Prevent index error
        % Calculate the time step (s)
        t_step = tout(t+1) - tout(t);
        pause(t_step);              % Assume negligible processing time
    end
   
end