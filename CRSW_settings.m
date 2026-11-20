function cfg = CRSW_settings()
% Configure one single-polarization data set.
cfg.input_dir = 'DGB';
cfg.output_dir = fullfile(pwd, 'CRSW_output');
cfg.master_date = '20170224';
cfg.nlines = 681;
cfg.slc_format = 'cpxint16';

cfg.crop = [];                 % [first row, last row, first column, last column]
cfg.patches = [1, 1];          % [azimuth, range]
cfg.window_size = 10;
cfg.overlap = 3;
cfg.shp_half_window = [5, 5];  % [azimuth, range]
cfg.shp_alpha = 0.05;
cfg.shp_prefilter = false;
cfg.min_shp = 1;
cfg.workers = 20;               % 0 runs the pixel loops serially
end
