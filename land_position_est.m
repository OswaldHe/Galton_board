function px = land_position_est(px0, R, py)
% estimate the landing position of a bead initially at x=x0
% the target y-displacement is assumed to be py
% The coordinates are relative to the center of the current peg
%
% p0:   initial x-coordinate
% R:    radius of peg
% py:   y-coordinate of the land position

%% input and constants
if px0 >= 0
    id = 1;
else
    px0 = -px0;
    id = 0;
end
g = 9.81;

%% theta_0 and v0
theta_0 = asin( px0/R );

% inelastic
v0 = 0;

%% calculate escaping velocity
[vt, theta_t] = output_velocity(v0, theta_0, R);

% initial position
x0 = R * sin(theta_t);
y0 = R * cos(theta_t);

% initial velocity
vx = vt * cos(theta_t);
vy = vt * sin(theta_t);


%% find the x-coordinate at y=py
% solve  1/2*g*t^2 + vy * t = py-y0
A = 1/2*g;
B = vy;
C = -(py-y0);
t = (-B + sqrt(B^2-4*A*C)) / (2*A);
px = x0 + vx*t;

if id == 0
    px = -px;
end

end
