%% Rahul_Goyal_main Usage and Description
% ME 326 Winter 2018 - Laboratory Assignment #4
%
% *Author:* RAHUL GOYAL
%
% California Polytechnic State University, San Luis Obispo, CA
%
% *Date Created:* February 06, 2018
%
% *Date Modified:* February 13, 2018
%
% *Description:*
% TODO
% 
% *Required Files:*
%
% * BowlingBallEOM.m - TODO
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



%% Initial Conditions
% TODO

% Initial Conditions
x_0 = [1.5;
       30;
       0;
       -25;
       0;
       0];



%% Simulink
% TODO

sim('example');



%% Plotting
% TODO

for i = 1
    
    % Plot the lane
    plot_lane();
    hold on;
    
    % TODO
    viscircles([0, 1], 1, 'Color', 'b')
    
end