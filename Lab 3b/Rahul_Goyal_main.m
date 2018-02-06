%% Rahul_Goyal_main Usage and Description
% ME 326 Winter 2018 - Laboratory Assignment #3B
%
% *Author:* RAHUL GOYAL
%
% California Polytechnic State University, San Luis Obispo, CA
%
% *Date Created:* January 30, 2018
%
% *Date Modified:* February 06, 2018
%
% *Description:*
% TODO
%
% *Required Files:*
%
% * CarEOM.m - TODO
% * TODO
%
% *Still To Do:*
%
% * Fix



%% Problem Statement
% TODO



%% Reset
% The following was used while debugging.

close all;
clear all;
clc;



%% Set Values
% The following is used to easily change the set horizontal velocity of the
% car. TODO?
v_x = 10;                       % Horizontal velocity of car (ft/s)



%% Given Values
% The following assigns values given by the problem statement to variables.
% (The value of the gravitational constant is assumed).

% Given Values
w_c = 750;                      % Weight of car's sprung mass (lb)
w_u = 85;                       % Weight of car's unsprung mass (lb)
k_s = 200;                      % Spring constant of the suspension (lb/in)
b_s = 60;                       % Damping coeff of the suspension (lb*s/in)
k_t = 1600;                     % Spring constant of tire (lb/in)
g = 32.174*12;                  % Gravitational constant (in/s^2)
h = 4;                          % Height of curb (in)
w = 6/12;                       % Width of curb (ft)

% Initial conditions (y_r, y_u, y_c, yu_dot, yc_dot)
ICs = [0;                       % Road displacement (in)
       0;                       % Car's unsprung mass dispacement (in)
       0;                       % Car's sprung mass displacement (in)
       0;                       % Car's unsprung mass velocity (in/s)
       0];                      % Car's sprung mass velocity (in/s)



%% Solved Values
% The following assigns values derived and/or solved from the given values
% to variables. See the attached sheet for hand calculations.

% Solved Values
m_c = w_c/g;                    % Mass of car's sprung mass (slinches)
m_u = w_u/g;                    % Mass of car's unsprung mass (slinches)

% State-Space Matrices
A = [0, 0, 0, 0, 0;
     0, 0, 0, 1, 0;
     0, 0, 0, 0, 1;
     k_t/m_u, -(k_s+k_t)/m_u, k_s/m_u, -b_s/m_u, b_s/m_u;
     0, k_s/m_c, -k_s/m_c, b_s/m_c, -b_s/m_c];

B = [1;
     0;
     0;
     0;
     0];

C = [1, 0, 0, 0, 0;
     0, 1, 0, 0, 0;
     0, 0, 1, 0, 0;
     0, 0, 0, 0, 0;
     0, 0, 0, 1, 0;
     0, 0, 0, 0, 1;
     0, -k_s, k_s, -b_s, b_s];

D = [0;
     0;
     0;
     1;
     0;
     0;
     0];



%% Generate a Curb Profile
% TODO COMMENT AND ORGANIZE SECTION

%%
% TODO EXPLAIN LARGE RESOLUTION = SMALLER STEPS
% TODO EXPLAIN TUNING: REALISM VS. PERFORMANCE
res = 5;                        % Resolution scalar (unitless)



%%
% TODO

% Curb Positions
d1 = 1;                         % Curb start position (ft)
d2 = d1+w;                      % Curb end position (ft)

% Slope Impulses
dw = 0.001;                     % Slope impulse "width" (ft)
k_w = h/dw;                     % Slope impulse "height" (in/ft)

% Curb Profile
d_f = 3/4*v_x;                  % Simulation distance (ft)
d = (0:dw/res:d_f)';            % Positions list (ft)
slope = [d, zeros(size(d))];	% Empty curb profile (in/ft)

%%
% The following four lines add a (positive) "impulse" at from d1 to d1+dx
% and a (negative) "impulse" from d2 to d2+dx to simulate a dirac delta
% function at the start and end conditions of the curb. The code does this
% by adding the calculated impulse for all values corresponding to greater
% than d1, then subtracting the calculated impulse for all values
% corresponding to greater than d1+dx. Thus, only the values between d1 and
% d1+dx have the calculated impulse added to them. A similar procedure, but
% negated, is followed to add the negative impulse to the values between d2
% and d2+dx.

% Add k_w to all corresponding slopes greater than d1
slope(d>=d1, 2) = slope(d>=d1, 2) + k_w;
% Subtract k_w from all corresponding slopes greater than d1 + dw
slope(d>=d1+dw, 2) = slope(d>=d1+dw, 2) - k_w;
% Add negative k_w to all corresponding slopes greater than d2
slope(d>=d2, 2) = slope(d>=d2, 2) - k_w;
% Subtract negative k_w from all corresponding slopes greater than d2 + dw
slope(d>=d2+dw, 2) = slope(d>=d2+dw, 2) + k_w;

% Debugging statement
% plot(slope);

%%
% The following TODO

t_step = (dw/v_x)/res;          % Time step (s)
t_f = d_f/v_x;                  % Time final (s)
t = (0:t_step:t_f)';            % Times list (s)



%% Simulating the Stuck Quarter-Car Using lsim
% The following TODO

% State-Space System Object
sys = ss(A, B, C, D);  

xr_stuck = v_x*t;               % Corresponding distance values (ft)
% Interpolate corresponding slope values from road_slope_data
mr_stuck = interp1(slope(:, 1), slope(:, 2), xr_stuck); % (in/ft)
u_stuck = v_x*mr_stuck;         % Corresponding u (yr_dot) values (in/s)

% Solve for y
y_stuck = lsim(sys, u_stuck, t);

% Debugging statement
% plot(xr_stuck, y_stuck);        % Plots y (7 variables) vs. x_r (stuck)



%% Simulating the Unstuck Quarter-Car Using ode45
% The following TODO

% TODO TEMP
options = odeset('MaxStep', dw/res);

% Quarter-Car Differential Equation Setup
CarODE = @(t, x) CarEOM(t, x, A, B, slope, v_x, k_t, m_u, w_c, w_u);
% Solve for x
[t, x_unstuck] = ode45(CarODE, t, ICs, options);

plot(t, x_unstuck(:, 1));
figure;
plot(t, x_unstuck(:, 2));
figure;
plot(t, x_unstuck(:, 3));
figure;
plot(t, x_unstuck(:, 4));
figure;
plot(t, x_unstuck(:, 5));




%% Discussion
% TODO
% TODO EXPLAIN LSIM VS. ODE45


