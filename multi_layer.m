%% demo for *****
% coordinate system
% centered at the center of the current circle
% x: to the right
% y: upwards

%% parameters
R = 1;
W = [1.881  1.881  3.12  3.117] * R;
H = 4 * R;
g = 9.81;
N_LV = length(W);

% % plot the two circles
figure;
hold on;
axis equal
theta_tmp = 0 : 359;
plot( cosd(theta_tmp), sind(theta_tmp) );
plot( cosd(theta_tmp)+W, sind(theta_tmp)-H );


%% random initial positions
% number of trials
num = 500;

% normally distributed initial positions
% normal_mean = 0;
% normal_std = pi/2 / 6;
% T_0 = random('Normal',normal_mean, normal_std, num, 1 );

% uniformly sampled
% x0 = R/2 /num : R/2 /num : R/2;     % single sided
x0 = -R/2 : R/(num-1) : R/2;     % double sided

% uniformly distributed initial positions
% x0 = random('Uniform', -R/2, R/2, 1, num);

% zero initial velocity
V0 = zeros(1,num);

% % non-zero initial velocity
% V0_r = 0.3;
% V0_max = sqrt( g * R * cos(T_0) );
% V0 = V0_max * V0_r;

%% calculate the landing positions
xt = zeros(1, num);     % landing positions
% yt = zeros(1, num);
pt_id = zeros(1, num);  % landing peg id

% initial
xi = x0;    % initial positions
p0_id = ones(1,num);	% peg id
p0 = 0;     % peg positions of current level

for lv = 1:N_LV
    % peg positions of the next level
    w = W(lv);
    if lv == 1
        pt = [-w  w];
    else
        pt = [ p0(1)-w  ...
               ( p0(2:end) + p0(1:end-1) )/2 ...
               p0(end)+w ];
    end    
    
    % landing positions of beads
    for i = 1 : num
        % target peg
        if xi(i) >= 0
            pt_id(i) = p0_id(i) + 1;
        else
            pt_id(i) = p0_id(i);
        end
        
        % relative distance of the current and the target peg
        w = pt( pt_id(i) ) - p0( p0_id(i) );
            
        % landing position (local coordinates)
        xt(i) = land_position( x0(i), R, [w H] );

        % landing position (global coordinates)
        xt(i) = xt(i) + p0( p0_id(i) );
        
        % update the initial position (prepare for next landing)
        xi(i) = xt(i) - pt( pt_id(i) );
    end
    
    % plot (for debugging)
    hist(xt);
   
    % update
    p0 = pt;
    p0_id = pt_id;

end

hist(xt);
% [sum(xt < W)  sum(xt>=W)]
% [ median(xt(xt>0))  median(xt(xt<0)) ]


