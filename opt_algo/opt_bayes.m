function [ opt_val, opt_pt, opt_hist, x_hist, calc_time_hist ] = opt_bayes( test_fn, x0, bounds, max_iter )
    if length( x0 ) == 5
        x1 = optimizableVariable('x1',bounds(1,:));
        x2 = optimizableVariable('x2',bounds(2,:));
        x3 = optimizableVariable('x3',bounds(3,:));
        x4 = optimizableVariable('x4',bounds(4,:));
        x5 = optimizableVariable('x5',bounds(5,:));
        vars = [ x1, x2, x3, x4, x5 ];
        init_x = array2table(normd2real( x0, bounds )');
    else
        x1 = optimizableVariable('x1',bounds(1,:));
        x2 = optimizableVariable('x2',bounds(2,:));
        x3 = optimizableVariable('x3',bounds(3,:));
        x4 = optimizableVariable('x4',bounds(4,:));
        x5 = optimizableVariable('x5',bounds(5,:));
        x6 = optimizableVariable('x6',bounds(6,:));
        x7 = optimizableVariable('x7',bounds(7,:));
        x8 = optimizableVariable('x8',bounds(8,:));
        x9 = optimizableVariable('x9',bounds(9,:));
        x10 = optimizableVariable('x10',bounds(10,:));
        vars = [ x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 ];
        init_x = array2table(normd2real( x0, bounds )');
    end
    res = bayesopt( test_fn, vars, 'MaxObjectiveEvaluations', max_iter, 'InitialX', init_x, 'Verbose', 0, 'PlotFcn', [], ...
        'AcquisitionFunctionName', 'expected-improvement', 'IsObjectiveDeterministic', true );
    opt_val = res.MinObjective;
    opt_pt = table2array(res.XAtMinObjective);
    opt_hist = res.ObjectiveMinimumTrace';
    x_hist = table2array(res.XTrace)';
    calc_time_hist = res.IterationTimeTrace';
    
