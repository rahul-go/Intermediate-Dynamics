%% boatODE Usage and Description
% Update this file with your own boat ODE. Remember to update the comments!
%
% The first function file contains the equations of motion in the form:
% function [xdot] = boatODE(t,x) Where t is time, x is the state vector and
% xdot is the derivative of the state vector. The only code contained in
% this function file should be the equations of motion in matrix form and
% whatever parameters are needed to define your EOM. Notice that the EOM
% Function, boatODE, utilizes the current state vector to evaluate the
% current time derivative of the state vector; this relationship will be
% key when solving state-determined problems. The above declaration is of
% the standard form to be solved by any of the built-in ODE solvers that
% MATLAB™ offers. That is, the function takes only two parameters t and x,
% in that order, and returns only xdot.
%
% Replace the following code with the equations of motion you determined
% during your hand calculations
function [xdot] = boatODE(t, x)

m = 147150 / 9.81;              % Fishing boat mass (kg)
k = 1000 / 0.0278;              % Spring constant (N/m)
b = 55000;                      % Drag coefficient (N*s/m)
v_t = 2;                        % Tug boat velocity (m/s)

x_t = v_t * t;                  % Tug boat position (x)

xdot = [ x(2)
         k/m*(x_t-x(1)) - b/m*x(2)];

end