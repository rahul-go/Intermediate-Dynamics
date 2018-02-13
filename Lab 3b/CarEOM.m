%% CarEOM Usage and Description
% This function represents the equations of motion for the simulation. It
% returns xdot with inputs of time, x, the A and B state-space matrices,
% the slope data, the horizontal velocity of the car, the tire's spring
% constant, the car's unsprung mass, the car's sprung weight, and the car's
% unsprung weight. See the attached sheet for hand calculations.

function [xdot] = CarEOM(t, x, A, B, slope, v_x, k_t, m_u, w_c, w_u)

%% Solve for u
% To find u (in this case, the vertical component of the velocity of the
% car) corresponding to the time and x value, the following solves for the
% corresponding distance value, from which it interpolates the
% corresponding slope value, from which it solves for the corresponding u
% value, the vertical component of the velocity of the car. Then, knowing
% the x and u values, the script solves y = A*x + B*u using matrix
% multiplication.

x_r = v_x*t;                    % Distance of car (ft)
% Interpolate slope (m) from distance (in/ft)
m_r = interp1(slope(:, 1), slope(:, 2), x_r);
u = v_x*m_r;                    % Vertical velocity (yr_dot) of car (ft/s)



%% Solve for xdot
% See the attached sheet for hand calculations.
xdot = A*x + B*u;



%% Unstick Tires
% To unstick the tires, the following checks whether the force in the
% tire's "spring" exceeds the maximum downward force before the tire
% separates from the ground. If it does, the acceleration of the unsprung
% mass must be updated. (The force does not act on the sprung mass). The
% maximum force is added, the calculated force is subtracted, and this sum
% is divided by mass to obtain acceleration. This acceleration is the
% "error" in the calculated acceleration of the unsprung mass. The error is
% added as a negated value to the acceleration of the unsprung mass because
% the error points in the opposite direction of the assumed direction of
% the acceleration of the unsprung mass.

y_r = x(1);                     % Easy access to y_r
y_u = x(2);                     % Easy access to y_u

f_kt = k_t*(y_u-y_r);           % Force in spring k_t (stuck) (N)
fkt_max = w_c + w_u;            % Saturation threshold (N)

if f_kt > fkt_max
    % Update the acceleration of the unsprung mass
    xdot(4) = xdot(4) - (fkt_max-f_kt)/m_u;
end
    
end