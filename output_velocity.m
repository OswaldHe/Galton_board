%% calculuate the output velocity
% v0:       initial velocity
% theta_0:  initial angle (relative to vertical axis);
% v_t:      output velocity
% theta_t:  output angle

% Updated on 2019/1/3


function [vt, theta_t] = output_velocity(v0, theta_0, R)

g = 9.81;

%% calculate the velocity
vt = sqrt( (v0^2 + 2 * g * R * cos(theta_0)) / 3 );  % output velocity

%% calculate the angle
theta_t = acos( (v0^2 + 2*g*R*cos(theta_0)) / (3*g*R) );       % output angle

