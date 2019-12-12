function [px,py] = land_position(px0, R, C)
% Calculate the landing position of a bead initially at x=x0
% The coordinates are relative to the center of the current peg
%
% p0:   initial x-coordinate
% R:    radius of peg
% C:    [W, H], coordinates of the center of the peg in the next level

%% input and constants
if px0 >= 0
    W = C(1);
    H = C(2);
    id = 1;
else
    px0 = -px0;
    W = -C(1);
    H = C(2);
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

% plot (for debugging)
% figure;
% hold on;
% axis equal
% theta_tmp = 0 : 359;
% plot( cosd(theta_tmp), sind(theta_tmp) );
% plot( cosd(theta_tmp)+W, sind(theta_tmp)-H );
% t_tmp = 0 : W/vx/50 : W/vx;
% plot( x0 + vx*t_tmp, y0-vy*t_tmp-1/2*g*t_tmp.^2 );

%% find the intersection with the next circle
% construct the polynomial of time t
poly_A = [ ( 1/4 * g^2 );
           ( vy * g );
           ( vx^2 + vy^2 - g*(y0+H) );
           ( 2 * vx * (x0-W) - 2 * vy * (y0+H) );
           ( (x0-W)^2 + (y0+H)^2 - R^2 )];

% find the roots of the polynomial
poly_A_roots = roots(poly_A);

% discard the complex roots
ind = abs( imag(poly_A_roots) ) < 1e-4;
poly_A_roots = poly_A_roots(ind);

if isempty(poly_A_roots)
    px = nan;
    py = nan;
else
    % find the travelled time
    t = min(poly_A_roots);

    % find the coordinates of the intersection
    px = x0 + vx * t;
    py = y0 - ( vy * t + 1/2 * g * t^2 );
    if id == 0
        px = -px;
    end
end



end
