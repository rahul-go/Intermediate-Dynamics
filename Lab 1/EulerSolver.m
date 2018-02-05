%% EulerSolver Usage and Description
% Update this file with your own boar ODE. Remember to update the comments!
%
% The second function file is a numerical integrator utilizing Euler’s
% method. The function declaration is of the form: function [t, x_array] =
% EulerSolver(func hand, t_f, x_0, t_step) Assignment 1-4 Where func_hand
% is a function handle for the ODE we are trying to integrate. In this case
% the handle should be @boatODE. This will give EulerSolver something to
% integrate. The final time, t_final, specifies when to end the simulation.
% All ODEs require initial conditions to solve; these will be passed in as
% x_0, a vector containing the values of the state variables at time zero.
% Finally, t_step determines the step size EulerSolver uses when computing
% the solution to the ODE. Refer to Appendix A on how to implement Euler’s
% method. The function EulerSolver should include creating a time array,
% defining initial conditions, and a loop implementing Euler’s method to
% successively integrate the ODE(s) passed in. Notice that Euler’s method
% is dependent on the current time derivative of the state vector to
% evaluate the state vector at the next point in time; an inverse
% relationship from the EOM function.
%
% Replace the following code with the EulerSolver algorithm discussed in
% the handout.

function [t, x_array] = EulerSolver(func_hand, t_f, x_0, t_step)

% Step 1: Create a column vector, t, containing the value of time for each
%  iteration of the solution, starting with 0, counting by t_step, and
%  ending at t_f.
t = (0:t_step:t_f)';            % Time

% Step 2: Create an array containing zeros, x_array, that is the same
%  number of rows as the time vector, t, and the same number of columns as
%  there are rows in initial state vector, x_0.
numt = length(t);               % Number of time points
numx = length(x_0);             % Number of variables
x_array = zeros(numt, numx);    % Initializing the x_array

% Step 3: Replace the first row of x_array with the initial conditions,
%  x_0, so that each following row can be computed.
x_array(1, :) = x_0;

% Step 4: Using a for loop, iterate through each row in x_array
for i = 1:numt-1

%  Step 4a: Compute the time rate of chanxge of the vector x by executing
%   the function_handle.
    xdot = func_hand(t(i), x_array(i, :));

%  Step 4b: Using the current time rate of change, the time step, and the
%   current value of x_array, compute the next value of x_array.
    x_array(i+1, :) = x_array(i, :) + xdot'*t_step;

end
end