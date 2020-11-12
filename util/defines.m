test_case     = 0;
FN_ROSEN      = 1;
FN_STYB_5D    = 2;
FN_STYB_10D   = 3;
FN_DEB1_5D    = 4;
FN_DEB1_10D   = 5;
FN_DEB2_5D    = 6;
FN_DEB2_10D   = 7;
FN_SCHWEF_5D  = 8;
FN_SCHWEF_10D = 9;
FN_SAL_5D     = 10;
FN_SAL_10D    = 11;
test_fn_names = { 'ROSEN', 'STYB_5D', 'STYB_10D', 'DEB1_5D', 'DEB1_10D', ...
                  'DEB2_5D', 'DEB2_10D', 'SCHWEF_5D', 'SCHWEF_10D', 'SAL_5D', 'SAL_10D' };

opt_tech      = 0;
OPT_RAND      = 1;
OPT_ADALIPO   = 2;
OPT_SMGO      = 3;
OPT_DIRECT    = 4;
OPT_BAYES     = 5;
opt_tech_names = { 'RAND', 'ADALIPO', 'SMGO', 'DIRECT', 'BAYES' };

FIG_OPTVALHIST = 1;
FIG_TESTPTHIST = 2;
FIG_OPTVALSTAT = 3;
FIG_VALHIST    = 4;
FIG_DEBUG1     = 5;
FIG_DEBUG2     = 6;

max_trial = 5;
max_iter = 10;

opt_params = { [ 0.0 ]; ...          % RAND
               [ 0.1 ]; ...          % ADALIPO
               [ 0.015 ]; ...        % SMGO (in this list are mu values, alpha is fixed at 0.015)
               [ 0.0 ]; ...          % DIRECT
               [ 0.0 ]; };           % BAYES