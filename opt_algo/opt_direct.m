function [ opt_val, opt_pt, opt_hist ] = opt_direct( test_fn, bounds, max_iter )
    prob.f = @(c) test_fn(c);
    opts.maxevals = max_iter;
    opts.maxits   = 1000;
%     opts.showits  = 0;
    [ opt_val, opt_pt, opt_hist ] = Direct( prob, bounds, opts );