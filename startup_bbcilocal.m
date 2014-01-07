global DATA_DIR
global TMP_DIR LOG_DIR

TMP_DIR= '/tmp/';
LOG_DIR= '/tmp/';
DATA_DIR= '/mnt/data/Master_Results/Datasets/bbciRaw/';

% avoid "assert overload" conflict:
warning('off','MATLAB:dispatcher:nameConflict');
% avoid annoying behaviour of set_defaults:
warning('off','set_defaults:DEFAULT_FLD');

startup_bbci;
