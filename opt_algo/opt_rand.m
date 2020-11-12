function [ opt_z, opt_pt, z_hist, opt_z_hist, x_hist, calc_time_hist ] = opt_rand( test_fn, x0_normd, bounds, max_iter )
% test_fn - test function
% x0 - starting point (N x 1) column matrix, can be []
% bounds - lower and upper bounds. must be an N x 2 matrix
% opt_param - bias parameter between exploration and exploitation
% max_iter - maximum iterations for a run

%% define initial variable values

opt_samp = [];
opt_z = Inf;

z_hist       = NaN * ones( 1, max_iter );
opt_z_hist   = NaN * ones( 1, max_iter );
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

    %% generate random point
    calc_time = tic;
    x_normd = rand( size(x0_normd) );
    calc_time_hist( iter ) = toc(calc_time);
end