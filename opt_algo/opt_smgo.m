function [ opt_z, opt_x, z_hist, opt_z_hist, x_hist, calc_time_hist ] = opt_smgo( test_fn, x0, bounds, alpha, mu, max_iter )
% test_fn - test function
% x0 - starting point (N x 1) column matrix, can be []
% bounds - lower and upper bounds. must be an N x 2 matrix
% opt_param - bias parameter between exploration and exploitation
% max_iter - maximum iterations for a run

%% Notes:
% - The algorithm should normalize the ranges to [0, 1] in order to
% generate a hypercube
[ D , ~ ] = size( bounds );

%% define initial variable values
gamma_0    = 1e-6;    
gamma      = gamma_0; 
gamma_prev = gamma_0; 

opt_samp = [];
opt_z    = Inf;

opt_z_hist     = NaN * ones( 1, max_iter );
calc_time_hist = NaN * ones( 1, max_iter );

% generate all corners of the hypercube (2^N)
vtx_hprbx = [];
for i = 0:2^D-1
    vtx_hprbx = [ vtx_hprbx dec2bin(i, D)'-'0' ];
end
vtx_corner_dist = D*ones( 1, 2^D );

%% Construct the database of segments and candidate points
% database of Mode theta: from all existing points (including corners),
% calculate the candidate point, and the exact lower bound at those points
% array structure:
% - rows 1 to D: the candidate point coordinate
% - row D+1: the value of (exact) lower bound cone vertex at that point
% - row D+2: the value of (exact) lower bound cone height at that point
% - row D+3: the value of (exact) lower bound at that point
THT_ROW_LBND_VTX = D + 1;
THT_ROW_LBND_HGT = D + 2;
THT_ROW_LBND     = D + 3;
mode_db_theta = [];

% database of Mode psi: from starting points, calculate the midpoint, and
% set the initial value of Lambda
% array structure:
% - rows 1 to D (end-5): the midpoint coordinate
% - row D+1 (end-4): the value of upper bound cone vertex at that midpoint
% - row D+2 (end-3): the value of upper bound cone height at that midpoint
% - row D+3 (end-2): the value of lower bound cone vertex at that midpoint
% - row D+4 (end-1): the value of lower bound cone height at that midpoint
% - row D+5 (end): the value of lambda at that midpoint
PSI_ROW_UBND_VTX = D + 1;
PSI_ROW_UBND_HGT = D + 2;
PSI_ROW_LBND_VTX = D + 3;
PSI_ROW_LBND_HGT = D + 4;
PSI_ROW_LAMBDA   = D + 5;
mode_db_psi = [];

%% loop to max iterations
for iter = 1:max_iter    
    if iter == 1
        opt_samp = [ vtx_hprbx; zeros( 1, 2^D ) ]; % initial values of the corners ( all with z = 0.0 )
        x = x0;        
    end
    
    z = test_fn( normd2real( x, bounds ) );
    opt_samp = [ opt_samp [ x; z ] ];          % stacking the column vector into collection of samples 
    
    %% mirroring scheme for corners
    vtx_corner_dist_tmp = sum( ( x*ones( 1, 2^D ) - opt_samp( 1:D, 1:2^D ) ) .^ 2 );
    vtx_corner_update   = vtx_corner_dist > vtx_corner_dist_tmp;
    vtx_corner_dist( :, vtx_corner_update ) = vtx_corner_dist_tmp( :, vtx_corner_update );
    opt_samp( end, vtx_corner_update ) = z;
    
    %% generate/update Lipschitz constant

    calc_time = tic;
    if iter > 1        
        % calculating Lipschitz constant
        gamma = max( gamma_prev, mu * max( abs( opt_samp( end, 1:end-1 ) - z ) ...
                  ./ sqrt( sum( ( opt_samp( 1:end-1, 1:end-1 )' - ones( 2^D+iter-1, 1 ) * x' )' .^ 2 ) ) ) );
    end
    
    %% updating the current best point
    if z < opt_z
        opt_z = z;
        opt_x = x;
        
        %% refresh mode_db_theta
        % recalculate locations of the candidate points
        idx_wo_opt     = opt_samp( end, (2^D+1):end ) > opt_z;
        idx_wo_opt_len = 2^D + sum( idx_wo_opt );
        idx_wo_opt     = [ logical( 1:2^D ) idx_wo_opt ]; % concatenating the corner indices (1 to 2^D)
        
        dx_vec         = ( opt_samp( 1:end-1, idx_wo_opt ) - opt_x * ones( 1, idx_wo_opt_len ) );
        dx_abs         = sqrt( sum( dx_vec .^ 2 ) ); % distances of opt point to all other sample points
        dz             = opt_samp( end, idx_wo_opt ) - opt_z;
        s_i            = dz ./ dx_abs;
        lb_dist_normd  = ( 1 - s_i / gamma ) / 2;
        mode_db_theta  = [ opt_x * ones( 1, idx_wo_opt_len ) + dx_vec .* ( ones( D, 1 ) * lb_dist_normd ); zeros( 3, idx_wo_opt_len ) ];
        
        % calculating the exact lower bounds
        for new_cdpt = 1:idx_wo_opt_len
            lbnd_lookup = [ opt_samp( end, 2^D + (1:iter) ); ...
                            gamma * sqrt( sum( ( opt_samp( 1:end-1, 2^D + (1:iter) ) - mode_db_theta( 1:D, new_cdpt ) * ones( 1, iter ) ) .^ 2 ) ) ];

            [ ~, lb_idx ] = max( lbnd_lookup( 1, : ) - lbnd_lookup( 2, : ) ); % calculating individual lower bound values        
            mode_db_theta( [ THT_ROW_LBND_VTX THT_ROW_LBND_HGT ], new_cdpt ) = [ lbnd_lookup( 1, lb_idx );
                                                                                 lbnd_lookup( 2, lb_idx ) ];
        end
            
    else
        %% update mode_db_theta
        
        % update existing entries
        mode_db_theta( THT_ROW_LBND_HGT, : ) = ( gamma / gamma_prev ) .* mode_db_theta( THT_ROW_LBND_HGT, : ); % rescaling due to gamma
        [ ~, old_cdpts ] = size( mode_db_theta );
        lbnd_lookup = [ z * ones( 1, old_cdpts ); ...
                        gamma * sqrt( sum( ( x * ones( 1, old_cdpts ) - mode_db_theta( 1:D, : ) ) .^ 2 ) ) ];
        lbnd_update = ( lbnd_lookup( 1, : ) - lbnd_lookup( 2, : ) ) > ...
                      ( mode_db_theta( THT_ROW_LBND_VTX, : ) - mode_db_theta( THT_ROW_LBND_HGT, : ));
        mode_db_theta( THT_ROW_LBND_VTX, lbnd_update ) = z;
        mode_db_theta( THT_ROW_LBND_HGT, lbnd_update ) = lbnd_lookup( end, lbnd_update );
        
        % add new candidate point
        s_i = ( z - opt_z ) / sqrt( sum( ( x - opt_x ) .^ 2 ) );        
        mode_db_theta_iter = [ ( opt_x + ( ( 1 - s_i / gamma ) / 2 ) * ( x - opt_x ) ); 0; 0; 0 ];
        lbnd_lookup = [ opt_samp( end, 2^D + (1:iter) ); ...
                        gamma * sqrt( sum( ( opt_samp( 1:end-1, 2^D + (1:iter) ) - mode_db_theta_iter( 1:D ) * ones( 1, iter ) ) .^ 2 ) ) ];
        [ ~, lb_idx ] = max( lbnd_lookup( 1, : ) - lbnd_lookup( 2, : ) ); % calculating individual lower bound values        
        mode_db_theta_iter( [ THT_ROW_LBND_VTX THT_ROW_LBND_HGT ] ) = [ lbnd_lookup( 1, lb_idx );
                                                                        lbnd_lookup( 2, lb_idx ) ];
        mode_db_theta = [ mode_db_theta mode_db_theta_iter ];
    end
    opt_z_hist( iter ) = opt_z;
    
    % filter all members whose lower bound cone vertex value are
    % coincident with the opt_z
    mode_db_theta_filter = ~( round( mode_db_theta( THT_ROW_LBND_VTX, : ), 4 ) > round( opt_z, 4 ) );
    mode_db_theta = mode_db_theta( :, mode_db_theta_filter );
    mode_db_theta( THT_ROW_LBND, : ) = mode_db_theta( THT_ROW_LBND_VTX, : ) - mode_db_theta( THT_ROW_LBND_HGT, : );
    
    %% iteratively update the midpoint database mode_db_psi
    % - should save the vertex value of the upper (lower) bound
    % - should also save the triangle value
    % - only the triangle value is scaled
    % - if the new incoming upper (lower) bound is tighter than existing
    %   - replace the vertex value
    %   - replace the triangle value
    if iter ~= 1
        mode_db_psi( [PSI_ROW_UBND_HGT PSI_ROW_LBND_HGT], : ) = ...
            ( gamma / gamma_prev ) .* mode_db_psi( [PSI_ROW_UBND_HGT PSI_ROW_LBND_HGT], : );
        [ ~ , mdpts ] = size( mode_db_psi );
        
        % calculating if i should update the upper (lower) bounds
        new_cone_height  = gamma * sqrt( sum( ( mode_db_psi( 1:D, : ) - x * ones( 1, mdpts ) ) .^ 2 ) );
        new_upper_bnd    = z + new_cone_height; 
        new_lower_bnd    = z - new_cone_height;
        upper_bnd_update = ( mode_db_psi( PSI_ROW_UBND_VTX, : ) + mode_db_psi( PSI_ROW_UBND_HGT, : ) ) > new_upper_bnd;
        lower_bnd_update = ( mode_db_psi( PSI_ROW_LBND_VTX, : ) - mode_db_psi( PSI_ROW_LBND_HGT, : ) ) < new_lower_bnd;
        
        % replacing the upper (lower) bounds data if needed
        mode_db_psi( PSI_ROW_UBND_VTX, upper_bnd_update ) = z;
        mode_db_psi( PSI_ROW_UBND_HGT, upper_bnd_update ) = new_cone_height( :, upper_bnd_update );
        mode_db_psi( PSI_ROW_LBND_VTX, lower_bnd_update ) = z;
        mode_db_psi( PSI_ROW_LBND_HGT, lower_bnd_update ) = new_cone_height( :, lower_bnd_update );
    end
    
    %% introduce additional midpoints to mode_db_psi
    mode_db_psi_iter = [ [ ( vtx_hprbx + x * ones( 1, 2^D ) ) / 2; zeros( 5, 2^D ) ] ...
                         [ ( opt_samp( 1:end-1, 2^D + (1:iter-1) ) + x * ones( 1, iter-1 ) ) / 2; zeros( 5, iter-1 ) ] ];
    [ ~, new_mdpts ] = size( mode_db_psi_iter );
    % calculating the upper (lower) bound information
    for new_mdpt = 1:new_mdpts
        ulbnd_lookup = [ opt_samp( end, 2^D + (1:iter) ); ...
                         gamma * sqrt( sum( ( opt_samp( 1:end-1, 2^D + (1:iter) ) - mode_db_psi_iter( 1:D, new_mdpt ) * ones( 1, iter ) ) .^ 2 ) ) ];
        [ ~, ub_idx ] = min( ulbnd_lookup( 1, : ) + ulbnd_lookup( 2, : ) ); % calculating individual upper bound values
        [ ~, lb_idx ] = max( ulbnd_lookup( 1, : ) - ulbnd_lookup( 2, : ) ); % calculating individual lower bound values        
        mode_db_psi_iter( PSI_ROW_UBND_VTX:PSI_ROW_LBND_HGT, new_mdpt ) = [ ulbnd_lookup( 1, ub_idx );
                                                                          ulbnd_lookup( 2, ub_idx );
                                                                          ulbnd_lookup( 1, lb_idx );
                                                                          ulbnd_lookup( 2, lb_idx ) ];
    end
    mode_db_psi = [ mode_db_psi mode_db_psi_iter ];
    mode_db_psi( PSI_ROW_LAMBDA, : ) = ( mode_db_psi( PSI_ROW_UBND_VTX, : ) + mode_db_psi( PSI_ROW_UBND_HGT, : ) ) ... % lambda is upper bound ...
                                     - ( mode_db_psi( PSI_ROW_LBND_VTX, : ) - mode_db_psi( PSI_ROW_LBND_HGT, : ) );    % minus the lower bound
    
    mode_db_psi                                 = unique( mode_db_psi', 'rows' )'; % there might be duplicate midpoints, so i get only the unique ones
    mode_db_psi( :,all(mode_db_psi(1:D,:)==x) ) = [];
    
    %% Algorithm (1) exploitation
    exploit_ok = 0;
    [ ~, exploit_idx ] = min( mode_db_theta( THT_ROW_LBND, : ) );
    if mode_db_theta( THT_ROW_LBND, exploit_idx ) < opt_z - alpha * gamma
        exploit_ok = 1;
        x = mode_db_theta( 1:D, exploit_idx );
    end
    % if Algorithm (1) fails
    if ~exploit_ok
        [ ~, explore_idx ]             = max( mode_db_psi( PSI_ROW_LAMBDA, : ) );
        x                              = mode_db_psi( 1:D, explore_idx );        
        mode_db_psi( : , explore_idx ) = [];
    end
    
    gamma_prev         = gamma;
    
    calc_time_hist( iter ) = toc( calc_time );
    z_hist = opt_samp( end, : );
    x_hist = opt_samp( 1:end-1, : );
end
