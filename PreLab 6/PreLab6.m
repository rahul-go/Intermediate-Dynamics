%% Reset
% The following was used while debugging.

close all;
clear all;
clc;



%% Script
% The following is the main script.
t2_0 = deg2rad(0);              % Angular position initial of link 2 (rad)
t3_0 = 0;                       % Angular position initial of link 3 (rad) [GUESS]
t4_0 = 0;                       % Angular position initial of link 4 (rad) [GUESS]
minimize = @(x) MyPosIC(t2_0, x);
x_0 = fminsearch(minimize, [t3_0, t4_0]);