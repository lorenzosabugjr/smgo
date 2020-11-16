clear all; close all;
addpath('test_fns', 'opt_algo', 'util');
if ~exist('out_dat', 'dir')
    mkdir('out_dat');
end
addpath('out_dat');
defines;

%% Automated run for all test functions, opt algos, and exploit weighting
for test_fn = [ FN_ROSEN FN_STYB_5D FN_STYB_10D FN_DEB1_10D FN_DEB2_10D FN_SCHWEF_10D FN_SAL_10D ]
    test_init;
    x0_list_normd = rand( opt_dim, max_trial ); % randomly generate all the starting points for all the trials
    for opt_tech = [ OPT_RAND OPT_DIRECT OPT_ADALIPO OPT_SMGO OPT_BAYES ]
        for opt_param = opt_params{ opt_tech }
            x_hist_stack         = [];
            val_hist_stack       = [];
            opt_val_hist_stack   = [];
            calc_time_hist_stack = [];

            fprintf("SOLVING %s: %s, parameter = %.4f\n", test_fn_names{test_fn}, opt_tech_names{opt_tech}, opt_param);
            if opt_tech == OPT_DIRECT
                max_trial_x = 1;
            else
                max_trial_x = max_trial;
            end
            
            for trial = 1:max_trial_x
                x0 = x0_list_normd( : , trial ); % x0 is a column vector

                if opt_tech == OPT_RAND
                    [ opt_val, opt_pt, val_hist, opt_val_hist, x_hist, calc_time_hist ] = ...
                        opt_rand( fn_in, x0, bounds, max_iter );
                elseif opt_tech == OPT_ADALIPO
                    [ opt_val, opt_pt, val_hist, opt_val_hist, x_hist, calc_time_hist ] = ...
                        opt_adalipo( fn_in, x0, bounds, opt_param, max_iter );                    
                elseif opt_tech == OPT_SMGO
                    alpha = opt_param;
                    mu    = 1.025;
                    [ opt_val, opt_pt, val_hist, opt_val_hist, x_hist, calc_time_hist ] = ...
                        opt_smgo( fn_in, x0, bounds, alpha, mu, max_iter );
                elseif opt_tech == OPT_DIRECT
                    [ opt_val, opt_pt, opt_val_hist ] = opt_direct( fn_in, bounds, max_iter );
                elseif opt_tech == OPT_BAYES
                    [ opt_val, opt_pt, opt_val_hist, x_hist ] = opt_bayes( fn_in, x0, bounds, max_iter );
                    calc_time_hist = [];
                end
                
                if ~( opt_tech == OPT_DIRECT )
                    x_hist_stack            = [ x_hist_stack; x_hist ];
                    if ~( opt_tech == OPT_BAYES )
                        val_hist_stack          = [ val_hist_stack; val_hist ];
                    end
                    opt_val_hist_stack      = [ opt_val_hist_stack; opt_val_hist ];
                    calc_time_hist_stack    = [ calc_time_hist_stack; calc_time_hist ];
                end
                fprintf("Trial %4d: %f\n", trial, opt_val);
            end

            if ~( opt_tech == OPT_DIRECT )
                % save important data about the optimization runs
                mat_name = sprintf("out_dat/%s-%s-%.4f.mat", opt_tech_names{opt_tech}, test_fn_names{test_fn}, opt_param);
                save( mat_name, 'x0_list_normd', 'opt_val_hist_stack', 'x_hist_stack', 'calc_time_hist_stack' );
            end
        end
    end
end

