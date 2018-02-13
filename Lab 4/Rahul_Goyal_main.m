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
% * Integrator.slx - TODO
% * BowlingBallEOM.m - TODO
% * plot_lane.m - TODO
%
% *Still To Do:*
%
% * Display x, y of no-slip
% * Smaller sampling rate
% * Solid bowling ball
% * Simulink end time
% * Better solution to chatter
% * Comment throughout



%% Problem Statement
% Consider a bowling ball thrown with an initial angular velocity w0 and an
% initial linear velocity v0 from the instant it makes contact with an
% alley lane. The velocity at the contact point, vC, is due to the velocity
% of the center of mass of the bowling ball, vB, and due to the relative
% velocity, vC/B, caused by the angular velocity of the bowling ball. The
% surface of the bowling alley will be defined as the xy-plane, where the
% y-axis extends down the alley and the x-axis extends toward the right
% gutter. The direction of the contact velocity is what determines the
% direction of the friction force. Angle theta is defined as the angle
% between the positive extension of the x-axis and the contact velocity
% vector. Therefore, when theta = pi/2 the ball will be traveling straight
% down the alley toward the pins.
% 
% The friction present between the ball and the smooth wooden lane is best
% approximated by a Coulomb friction model. That is |F| = mukN if the ball
% is slipping and |F| <= musN if the ball is rolling without slip. In each
% regime the friction force opposes the relative velocity (or impending
% relative velocity) at the contact point. In order to develop equations of
% motion, the magnitude and direction of the coulomb friction force must be
% computed, and then put into component form in the xy plane.
% 
% The bowling ball is modeled as a sphere of uniform density that weighs 15
% lb and has a diameter of 8.5 in. The kinetic friction coefficient, muk,
% between the ball and surface is assumed to be 0.12 and the static
% friction coefficient, mus, between the ball and surface is assumed to be
% 0.14.



%% Reset
% The following was used while debugging.

close all;
clear all;
clc;



%% Initial Conditions
% TODO

% Initial Conditions
x_0 = [1.5;                     % Velocity[x] of ball (ft/s)
       30;                      % Velocity[y] of ball (ft/s)
       0;                       % Angular velocity[x] of ball (rad/s)
       -25;                     % Angular velocity[y] of ball (rad/s)
       0;                       % Displacement[x] of ball (ft)
       0];                      % Displacement[y] of ball (ft)



%% Simulate the Bowling Ball Using Simulink
% TODO

sim('Integrator');



%% Animation
%  TODO

for t = 1:length(tout)
    
    % If the bowling ball has traveled the lane length, finish
    lane_length = 62+(10+3/16)/12;  % Lane length (ft)
    if xout(t, 6) > lane_length
        break;
    end
    
    
    
    % Plot the Lane
    plot_lane();                    % TODO
    ylim([0, lane_length]);         % TODO
    axis equal;                     % TODO
    
    % Plot the Bowling Ball
    d = 8.5;                        % Given diameter (in)
    r = (d/2)/12;                   % Solved radius (ft)
    x = xout(t, 5);                 % X position of ball's center (ft)
    y = xout(t, 6);                 % Y position of ball's center (ft)
    
    % Acceleration[linear, angular] (ft/s, ft/s, rad/s^2, rad/s^2)
    a = xdot(t, 1:4);
    
    % If no slip, set bowling ball color to green
    if a == 0
        color = 'green';
    % Else, set bowling ball color to blue
    else
        color = 'blue';
    end
    
    % Draw circle (bowling ball)
    viscircles([x, y], r, 'Color', color);
    
    
    
    % Calculate the time step and pause accordingly
    if t ~= length(tout)            % Prevent index error
        % Calculate the time step (s) and store for later use
        t_step(t) = tout(t+1) - tout(t);
        pause(t_step(t));           % Assume negligible processing time
    end
    
    
    
    % Format Plot
%     title('Bowling Ball Animation');
%     xlabel({'X Position (ft)'
%             ''
%             % Figure label
%             '\bfFigure 1: \rmBowling Ball Animation'});
%     ylabel('Y Position (ft)');
    
    
    
    % Convert plot frame to image and store for later use
    image{t} = frame2im(getframe());
    
end



%% Export as GIF
% TODO
% COMMENT ALL

file_name = 'BowlingBallAnimation.gif';
for i = 1:length(image)
    % TODO
    [A, map] = rgb2ind(image{i}, 256);
    % If first iteration, run setup code
    if i == 1
        imwrite(A, map, file_name, ...
                'LoopCount', inf, ...
                'DelayTime', t_step(i));
    % Else, append images
    else
        imwrite(A, map, file_name, ...
                'WriteMode', 'append', ...
                'DelayTime', t_step(i));
    end
end