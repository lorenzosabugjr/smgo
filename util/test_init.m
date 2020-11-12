if test_fn == FN_STYB_5D
    opt_dim = 5;
    bounds = ones( opt_dim, 1 ) * [ -5 5 ] ;
    fn_in = @( c ) fn_styblinski( c );
elseif test_fn == FN_STYB_10D
    opt_dim = 10;
    bounds = ones( opt_dim, 1 ) * [ -5 5 ] ;
    fn_in = @( c ) fn_styblinski( c );
elseif test_fn == FN_ROSEN
    opt_dim = 10;
    bounds = ones( opt_dim, 1 ) * [ -40 5 ] ;
    fn_in = @( c ) fn_rosenbrock( c );
elseif test_fn == FN_DEB1_5D
    opt_dim = 5;
    bounds = ones( opt_dim, 1 ) * [ -1 1 ] ;
    fn_in = @( c ) fn_deb1( c );
elseif test_fn == FN_DEB1_10D
    opt_dim = 10;
    bounds = ones( opt_dim, 1 ) * [ -1 1 ] ;
    fn_in = @( c ) fn_deb1( c );
elseif test_fn == FN_SCHWEF_5D
    opt_dim = 5;
    bounds = ones( opt_dim, 1 ) * [ -500 500 ] ;
    fn_in = @( c ) fn_schwefel( c );
elseif test_fn == FN_SCHWEF_10D
    opt_dim = 10;
    bounds = ones( opt_dim, 1 ) * [ -500 500 ] ;
    fn_in = @( c ) fn_schwefel( c );
elseif test_fn == FN_DEB2_5D
    opt_dim = 5;
    bounds = ones( opt_dim, 1 ) * [ 0 150 ] ;
    fn_in = @( c ) fn_deb2( c );
elseif test_fn == FN_DEB2_10D
    opt_dim = 10;
    bounds = ones( opt_dim, 1 ) * [ 0 150 ];
    fn_in = @( c ) fn_deb2( c );
elseif test_fn == FN_SAL_5D
    opt_dim = 5;
    bounds = ones( opt_dim, 1 ) * [ -40 70 ] ;
    fn_in = @( c ) fn_salomon( c );
elseif test_fn == FN_SAL_10D
    opt_dim = 10;
    bounds = ones( opt_dim, 1 ) * [ -40 70 ] ;
    fn_in = @( c ) fn_salomon( c );
end
