%% Reset
% The following was used while debugging.

close all;
clear all;
clc;



%% Script
% The following is the main script.

t2_0 = 0;                       % Theta 2 initial (rad)

t3_0 = 0;                       % Theta 3 initial (rad) [GUESS]
t4_0 = 0;                       % Theta 4 initial (rad) [GUESS]
x_0 = [t3_0, t4_0];

minimize = @(test_x) MyPosIC(t2_0, test_x);
x = fminsearch(minimize, x_0);