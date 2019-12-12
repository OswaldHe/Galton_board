%% beads released one after another into a Galton board
% The collision between the bead and the peg is assumed to be perfectly
% inelastic, i.e., the initial velocity is always zero
%
% The results show that
% to obtain a binomial distribution at the final level, the pegs must be
% precisely placed at particular positions.
% More concretely, a peg in the next level must be allocated where the
% beads landing on it are divided into two equal groups.
%
% It is also shown that after a few levels, the landing positions on a same
% peg tend to coincide, which makes it almost impossible to divide into
% equal groups.
%
% Yihao Zhou
% Last updated March 4, 2019


%% parameters
% radius of peg
R = 1;

% vertical distance between successive levels
H1 = 5 * R;     % between first two level
H = 4 * R;      % subsequent levels

% gravitational acceleration
g = 9.81;

% number of levels
N_LV = 6;

% to store peg positions in each level
PT = cell(1, N_LV);

%% random initial positions
% number of trials
num = 8096;

% normally distributed initial positions
% normal_mean = 0;
% normal_std = R/2 / 3;
% x_init = random('Normal',normal_mean, normal_std, [1 num] );
% x_init( x_init>R/2 ) = R/2;
% x_init( x_init<-R/2 ) = -R/2;

% uniformly sampled
x_init = -R/2 : R/(num-1) : R/2;     % double sided

%% calculate the landing positions
xt = zeros(1, num);     % landing positions
yt = zeros(1, num);
pt_id = zeros(1, num);  % landing peg id

% plot for debugging
figure;
    
for lv = 1:N_LV-1
    %% bead and peg positions in current level
    if lv == 1
        p0 = 0;     % peg at the center
        p0_id = ones(1, num);   % all beads fall on the only peg
        x0 = x_init;
        PT{1} = p0;
        h = H1;
    else
        p0 = pt;
        p0_id = pt_id;
        % local coordinates in current level
        x0 = xt - p0( p0_id );
        h = H;
    end
    
    %% peg positions in next level
    pt = zeros(1, lv+1);
    
    % interior pegs
    for k = 2 : lv
        % peg to the left in previous level
        id1 = find( p0_id == k-1  &  x0 > 0 );
        n1 = length(id1);
        
        % peg to the right in previous level
        id2 = find( p0_id == k  &  x0 < 0 );
        n2 = length(id2);
        
        if n1/n2 > 0.95  &&  n1/n2 < 1.05
            % approximately same
            pt(k) = ( PT{lv}(k-1) + PT{lv}(k) )/2;
        elseif n1 > n2
            % find the bead landing in the middle
            [~,sid] = sort(x0(id1), 'descend');
            mid = ceil((n1+n2)/2);
            id = id1( sid(mid) );
            % initial estimate
            pt(k) = land_position_est(x0(id), R, h) + PT{lv}(k-1);
            % refine
            iter = 1;
            while iter <= 3
                w = pt(k) - PT{lv}(k-1);
                pt(k) = land_position(x0(id), R, [w h]) + PT{lv}(k-1);
                iter = iter + 1;
            end
        else
            % find the bead landing in the middle
            [~,sid] = sort( x0(id2), 'ascend' );
            mid = ceil((n1+n2)/2);
            id = id2( sid(mid) );
            % initial estimate
            pt(k) = land_position_est(x0(id), R, h) + PT{lv}(k);
            % refine
            iter = 1;
            while iter <= 3
                w = pt(k) - PT{lv}(k);
                pt(k) = land_position(x0(id), R, [w h]) + PT{lv}(k);
                iter = iter + 1;
            end
        end
    end
    
    % rightmost peg
    id = find(p0_id==lv & x0>0);
    [~,sid] = sort( x0(id) );
    id = id( sid(floor(end/2)) );
    % initial estimate
    w = land_position_est(x0(id), R, h);
    % refine
    iter = 1;
    while iter <= 3
        w = land_position(x0(id), R, [w h]);
        iter = iter + 1;
    end
    pt(lv+1) = PT{lv}(lv) + w;
    
    % leftmost peg (symmetric)
    pt(1) = - pt(lv+1);
    
    % save
    PT{lv+1} = pt;
    
    %% landing positions
    for i = 1 : num
        % target peg
        if x0(i) >= 0
            pt_id(i) = p0_id(i) + 1;
        else
            pt_id(i) = p0_id(i);
        end
        
        % relative distance of the current and the target peg
        w = pt( pt_id(i) ) - p0( p0_id(i) );
            
        % landing position (local coordinates)
        [xt(i), yt(i)] = land_position( x0(i), R, [w h] );

        % landing position (global coordinates)
        xt(i) = xt(i) + p0( p0_id(i) );
        yt(i) = yt(i);
    end

%     % random noise
%     xt = xt + random('Normal',0,0.000001,[1 num]);
    
    % plot (for debugging)
%     hold off;
%     plot( xt, yt, 'r.','MarkerSize',12 );
%     hold on
%     axis equal
%     theta_tmp = 0 : 359;
%     for i = 1:lv+1 
%         plot( cosd(theta_tmp)+pt(i), sind(theta_tmp)-h );
%         plot( [pt(i) pt(i)],[-h+1.5*R, -h-1.5*R],'k' )
%     end
%     pause
    
    %% classify into bins
    bins = zeros(1,lv+1);     % number of beads in the bins

    % leftmost bin
    id = xt <= p0(1);
    bins(1) = sum(id);

    % interior bins
    for k = 2 : lv
        lower_xt = p0(k-1);
        upper_xt = p0(k);
        id = ( xt>lower_xt  &  xt<=upper_xt );
        bins(k) = sum(id);
    end

    % rightmost bin   
    id = xt > p0(lv);
    bins(lv+1) = sum(id);

    % plot (for debugging)
    fprintf('Level %d:',lv);
    fprintf('\t%.2f',bins/bins(1));
    fprintf('\n');
%     hist(bins/bins(1),lv+1);
    
   
end



