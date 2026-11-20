function manifest = CRSW_manifest(cfg)
% Match all input files by acquisition date.
required = {'input_dir','output_dir','master_date','nlines','slc_format', ...
    'crop','patches','window_size','overlap','shp_half_window', ...
    'shp_alpha','shp_prefilter','min_shp','workers'};
for k = 1:numel(required)
    assert(isfield(cfg, required{k}), 'CRSW:Settings', ...
        'Missing setting: %s.', required{k});
end
cfg.input_dir = char(cfg.input_dir);
cfg.output_dir = char(cfg.output_dir);
master = char(string(cfg.master_date));
assert(~isempty(regexp(master, '^\d{8}$', 'once')), 'CRSW:Date', ...
    'master_date must be a date in YYYYMMDD format.');
try
    parsed = datetime(master, 'InputFormat', 'yyyyMMdd', 'Format', 'yyyyMMdd');
    assert(strcmp(char(parsed), master));
catch
    error('CRSW:Date', 'Invalid master_date: %s.', master);
end
validateattributes(cfg.nlines, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.window_size, {'numeric'}, {'scalar','integer','>=',4});
validateattributes(cfg.overlap, {'numeric'}, ...
    {'scalar','integer','>=',1,'<',cfg.window_size});
validateattributes(cfg.shp_half_window, {'numeric'}, ...
    {'vector','numel',2,'integer','positive'});
validateattributes(cfg.patches, {'numeric'}, {'vector','numel',2,'integer','positive'});
validateattributes(cfg.shp_alpha, {'numeric'}, {'scalar','real','>',0,'<',1});
validateattributes(cfg.min_shp, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.workers, {'numeric'}, {'scalar','integer','nonnegative'});
assert(isscalar(cfg.shp_prefilter) && ...
    (islogical(cfg.shp_prefilter) || isnumeric(cfg.shp_prefilter)) && ...
    any(cfg.shp_prefilter == [0,1]), 'CRSW:Settings', ...
    'shp_prefilter must be true or false.');
assert(cfg.min_shp <= prod(2 * cfg.shp_half_window + 1), ...
    'CRSW:Settings', 'min_shp exceeds the number of pixels in the window.');
format = lower(char(cfg.slc_format));
assert(any(strcmp(format, {'cpxint16','cpxfloat32'})), 'CRSW:Settings', ...
    'slc_format must be cpxint16 or cpxfloat32.');
assert(isfolder(cfg.input_dir), 'CRSW:Input', 'Input directory not found: %s', cfg.input_dir);
assert(~isempty(cfg.output_dir), 'CRSW:Settings', 'output_dir must not be empty.');

flat_intensity = dir(fullfile(cfg.input_dir, '*.rmli'));
nested_intensity = dir(fullfile(cfg.input_dir, 'rmli', '*.rmli'));
assert(isempty(flat_intensity) || isempty(nested_intensity), 'CRSW:Input', ...
    'Intensity files exist in both the input directory and its rmli subdirectory.');
if ~isempty(flat_intensity)
    intensity = flat_intensity;
    slc_dir = cfg.input_dir;
    diff_dir = cfg.input_dir;
else
    intensity = nested_intensity;
    slc_dir = fullfile(cfg.input_dir, 'slc');
    diff_dir = fullfile(cfg.input_dir, 'diff');
end
assert(~isempty(intensity), 'CRSW:Input', ...
    'No .rmli files found in the input directory or its rmli subdirectory.');
dates = cell(1, numel(intensity));
for k = 1:numel(intensity)
    token = regexp(intensity(k).name, '^(\d{8})\.rmli$', 'tokens', 'once');
    assert(~isempty(token), 'CRSW:Input', 'Invalid intensity filename: %s', intensity(k).name);
    dates{k} = token{1};
    try
        date_value = datetime(dates{k}, 'InputFormat', 'yyyyMMdd', 'Format', 'yyyyMMdd');
        assert(strcmp(char(date_value), dates{k}));
    catch
        error('CRSW:Date', 'Invalid acquisition date: %s.', dates{k});
    end
end
[dates, order] = sort(dates);
intensity = intensity(order);
assert(numel(unique(dates)) == numel(dates), 'CRSW:Input', 'Duplicate acquisition dates.');
assert(numel(dates) >= 4, 'CRSW:Input', 'CCMmtf requires at least four acquisitions.');
master_index = find(strcmp(dates, master));
assert(numel(master_index) == 1, 'CRSW:Input', 'Master date %s is absent.', master);

width = intensity(1).bytes / (4 * cfg.nlines);
assert(width >= 1 && width == fix(width), 'CRSW:Input', ...
    'Intensity file size is inconsistent with nlines.');
slc_bytes = 4;
if strcmp(format, 'cpxfloat32'), slc_bytes = 8; end
slc = cell(size(dates));
diff = cell(size(dates));
for k = 1:numel(dates)
    assert(intensity(k).bytes == 4 * cfg.nlines * width, 'CRSW:Input', ...
        'Unexpected image size: %s.', intensity(k).name);
    candidates = {fullfile(slc_dir, [dates{k}, '.rslc']), ...
                  fullfile(slc_dir, [dates{k}, '.rslc_f'])};
    present = cellfun(@isfile, candidates);
    assert(sum(present) == 1, 'CRSW:Input', ...
        'Expected one SLC for %s (.rslc or .rslc_f).', dates{k});
    slc{k} = candidates{present};
    fileinfo = dir(slc{k});
    assert(fileinfo.bytes == slc_bytes * cfg.nlines * width, 'CRSW:Input', ...
        'Unexpected image size: %s.', slc{k});

    candidates = {fullfile(diff_dir, [dates{k}, '.diff']), ...
                  fullfile(diff_dir, [master, '_', dates{k}, '.diff'])};
    present = cellfun(@isfile, candidates);
    assert(sum(present) <= 1 && (k == master_index || any(present)), 'CRSW:Input', ...
        'Expected one master-to-date interferogram for %s.', dates{k});
    if any(present)
        diff{k} = candidates{present};
        fileinfo = dir(diff{k});
        assert(fileinfo.bytes == 8 * cfg.nlines * width, 'CRSW:Input', ...
            'Unexpected image size: %s.', diff{k});
    else
        diff{k} = '';
    end
end
if isempty(cfg.crop)
    crop = [1, cfg.nlines, 1, width];
else
    validateattributes(cfg.crop, {'numeric'}, {'vector','numel',4,'integer','positive'});
    crop = reshape(cfg.crop, 1, []);
    assert(crop(1) <= crop(2) && crop(3) <= crop(4) && ...
        crop(2) <= cfg.nlines && crop(4) <= width, 'CRSW:Settings', 'Invalid crop bounds.');
end
assert(cfg.patches(1) <= crop(2)-crop(1)+1 && cfg.patches(2) <= crop(4)-crop(3)+1, ...
    'CRSW:Settings', 'The number of patches exceeds the image dimensions.');
manifest = struct('dates', {dates}, 'master_date', master, 'master_index', master_index, ...
    'slc', {slc}, 'diff', {diff}, 'width', width, 'crop', crop);
end
