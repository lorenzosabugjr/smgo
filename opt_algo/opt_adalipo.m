function [ opt_z, opt_pt, z_hist, opt_z_hist, x_hist, calc_time_hist ] = opt_adalipo( test_fn, x0_normd, bounds, opt_param, max_iter )
% test_fn - test function
% x0 - starting point (N x 1) column matrix, can be []
% bounds - lower and upper bounds. must be an N x 2 matrix
% opt_param - bias parameter between exploration and exploitation
% max_iter - maximum iterations for a run

%% Notes:
% - The algorithm should normalize the ranges to [0, 1] in order to
% generate a hypercube
[ opt_dim , ~ ] = size( bounds );
bounds_normd = ones( opt_dim, 1 ) * [ 0 1 ];

%% define initial variable values
cone_slope_0 = 1e-6;    cone_slope   = cone_slope_0; 
cone_tol     = 0;       cone_defined = 0;

opt_samp = [];
opt_z = Inf;

z_hist         = NaN * ones( 1, max_iter );
opt_z_hist     = NaN * ones( 1, max_iter );
x_hist         = NaN * ones( length( x0_normd ), max_iter );
calc_time_hist = NaN * ones( 1, max_iter );

%% loop to max iterations
for iter = 1:max_iter
    
    if iter == 1
        x_normd = x0_normd; % TODO: must check if x0 is within bounds
    end
    x_hist( : , iter ) = x_normd;
    
    z = test_fn( normd2real( x_normd, bounds ) );
    opt_samp = [ opt_samp [ x_normd; z ] ]; % stacking the column vector into collection of samples
    
    if z < opt_z
        opt_z = z;
        opt_pt  = x_normd;
    end
    z_hist( iter ) = z;
    opt_z_hist( iter ) = opt_z;
    %% calculating next point to explore (AdaLIPO method)


    %% generate max slope

    calc_time = tic;
    data_length = length( opt_samp( 1, : ) );

    if data_length > 1
        cone_slope = max( cone_slope, max( abs( opt_samp( end, 1:end-1 ) - opt_samp( end, end ) ) ...
            ./  sqrt( sum( ( opt_samp( 1:end-1, 1:end-1 )' - ones( data_length-1, 1 ) * x_normd' )' .^ 2 ) ) ) );
    end
    log_1_5_slope = log(cone_slope)/log(1.5);
    cone_slope_adj = 1.5^( ceil( log_1_5_slope ) );

    %% generate cones and perform randomized selection (exploration-exploitation)

    if binornd( 1, opt_param )
        % exploit code here
        exploit_valid_flag = 0;
        exploit_trials = 0;
        while (exploit_valid_flag == 0) && (exploit_trials < 25000)
            x_normd = rand( size(x0_normd) );
            if lower_bnd( x_normd, opt_samp, cone_tol, cone_slope_adj ) < opt_z
                exploit_valid_flag = 1;
            end
            exploit_trials = exploit_trials + 1;
        end    
    else
        % explore code here
        x_normd = rand( size(x0_normd) );
    end

    calc_time_hist( iter ) = toc(calc_time);
end