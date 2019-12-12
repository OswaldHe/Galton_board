%% demo for *****
% coordinate system
% centered at the center of the current circle
% x: to the right
% y: upwards

%% parameters
R = 1;
D = 2.1 * R;
H = 5 * R;
g = 9.81;

% plot the two circles
% figure;
% hold on;
% axis equal
% theta_tmp = 0 : 359;
% plot( cosd(theta_tmp), sind(theta_tmp) );
% plot( cosd(theta_tmp)+D, sind(theta_tmp)-H );
% plot( cosd(theta_tmp)-D, sind(theta_tmp)-H );

%% generate theta_0
% number of trials
num = 300;
level = 4;
xt = 0;
yt = 0;
% normally distributed initial positions
% normal_mean = 0;
% normal_std = pi/2 / 6;
% T_0 = random('Normal',normal_mean, normal_std, num, 1 );

% uniformly distributed initial positions
T_0 = random('Uniform', -pi/3, pi/3, 1, num);

T_t = zeros(1, num);
for i = 1 : num
    %% theta_0
    theta_0 = T_0(i);
    for j = 1 : level
        
        if theta_0 >= 0
            flag = 1;
        else
            theta_0 = abs(theta_0);
            flag = 0;
        end

        %% calculate output velocity
        [v0, theta_c] = output_velocity(theta_0, R);

        % initial position
        %if flag == 1
            x0 = R * sin(theta_c);
            y0 = R * cos(theta_c);

            % initial velocity
            vx = v0 * cos(theta_c);
            vy = v0 * sin(theta_c);
%         else
%             x0 = -R * sin(theta_c);
%             y0 = R * cos(theta_c);
% 
%             % initial velocity
%             vx = -v0 * cos(theta_c);
%             vy = v0 * sin(theta_c);
%         end
        %% plot
%         t_tmp = 0 : H/vy/50 : H/3/vy;
%         plot( x0 + vx*t_tmp, y0-vy*t_tmp-1/2*g*t_tmp.^2 );
%         plot( x0 - vx*t_tmp, y0-vy*t_tmp-1/2*g*t_tmp.^2 );
%         plot( xt + vx*t_tmp, yt-vy*t_tmp-1/2*g*t_tmp.^2 );
        %plot( xt - vx*t_tmp, yt-vy*t_tmp-1/2*g*t_tmp.^2 );
    
        %% find the intersection with the next circle
        % construct the polynomial of time t
        poly_A = [ ( 1/4 * g^2 );
                   ( vy * g );
                   ( vx^2 + vy^2 - g*(y0+H) );
                   ( 2 * vx * (x0-D) - 2 * vy * (y0+H) );
                   ( (x0-D)^2 + (y0+H)^2 - R^2 )];

        % find the roots of the polynomial
        poly_A_roots = roots(poly_A);

        % discard the complex roots
        ind = abs( imag(poly_A_roots) ) < 1e-4;
        poly_A_roots = poly_A_roots(ind);

        if isempty(poly_A_roots)
            continue;
        end

        % find the time
        t = min(poly_A_roots);

        % find the coordinates of the intersection
        xt = x0 + vx * t;
        yt = y0 - ( vy * t + 1/2 * g * t^2 );
        
        %% save result and calculate new theta_0
        
        x0 = R * sin(theta_0);
        y0 = R * cos(theta_0);
        if flag == 1
            
            theta_0 = atan((xt-D)/(yt+H));
            if j==1
                T_t(i) = T_t(i) + xt;
            else
                T_t(i) = T_t(i) + xt-x0;
            end
        else
            theta_0 = -atan((xt-D)/(yt+H));
            if j==1
                T_t(i) = T_t(i) - xt;
            else
                T_t(i) = T_t(i) -(xt-x0);
            end
        end
    end
     histogram(T_t,5);
end


