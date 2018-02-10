function [ ] = plot_lane( )
%PLOT_LANE Plot the outline of a bowling lane with pins
%    PLOT_LANE() plots the shape of a bowling lane with a wood colored
%    background and all ten bowling pins as seen from the top. The pins
%    appear as circles filled in with white. The dimensions of the lane
%    conform to the standard for bowling alleys. That is the lane has
%    dimensions 42" W by 62'10-3/16"L. The bowling pins have a maximum
%    diameter of 4.766" as is standard. The pins are arranged in a
%    triangular grid, with 12" between centers. The center of the head-pin
%    is 60' from the beginning of the lane.

% For speed the constants necessary to plot the lane and pins are stored as
% persistent data.
persistent lane_length single_pin pin_array lane_array

% Generate the persistent data only if the data has not already been
% created
if isempty(pin_array)
    % Lane length in units of feet.
    lane_length = 62+(10+3/16)/12;
    
    % Array of data points forming a circle to represent a single pin.
    % Units are also in feet.
    single_pin = 4.766/2/12*[cos(0:pi/50:2*pi)' sin(0:pi/50:2*pi)'];
    
    % Cell array containing all 10 pins translated to their location
    % relative to the head-pin. This allows them each to be plotted
    % individually.
    pin_array = {single_pin+[ 0.0 0.0        ];   %Pin 1
                 single_pin+[-0.5 sqrt(3)/2  ];   %Pin 2
                 single_pin+[ 0.5 sqrt(3)/2  ];   %Pin 3
                 single_pin+[-1.0 sqrt(3)    ];   %Pin 4
                 single_pin+[ 0.0 sqrt(3)    ];   %Pin 5
                 single_pin+[ 1.0 sqrt(3)    ];   %Pin 6
                 single_pin+[-1.5 3*sqrt(3)/2];   %Pin 7
                 single_pin+[-0.5 3*sqrt(3)/2];   %Pin 8
                 single_pin+[+0.5 3*sqrt(3)/2];   %Pin 9
                 single_pin+[+1.5 3*sqrt(3)/2];}; %Pin 10
    
    % Array of points forming the outline of the lane.
    lane_array = [-21/12 0
                  -21/12 lane_length
                   21/12 lane_length
                   21/12 0
                  -21/12 0 ];
end

% This variable determines whether or not hold should be turned back on
% after the lane is done plotting.
hold_state = ishold;

% Plot the bowling lane itself and set the color to appear as wood.
fill(lane_array(:,1),lane_array(:,2),[240/255 190/255 150/255]);

% Plot the pins one by one as white circles translated to their location on
% the lane.
hold on;
for n=1:length(pin_array)
    fill(pin_array{n}(:,1),pin_array{n}(:,2)+60,'white');
end

% If needed, re-enable the plot hold so that further data (like the ball)
% may be plotted.
if ~hold_state
    hold off;
end
          
end

