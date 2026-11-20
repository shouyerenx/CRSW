function data = CRSW_read_SARimage(filename, nlines, width, format, bounds)
% Read a row-major big-endian raster block.
assert(bounds(1) >= 1 && bounds(2) <= nlines && bounds(1) <= bounds(2) && ...
    bounds(3) >= 1 && bounds(4) <= width && bounds(3) <= bounds(4), ...
    'CRSW:Input', 'Raster block is outside the image.');
switch lower(format)
    case 'cpxint16'
        precision = 'int16'; bytes = 2; components = 2;
    case 'cpxfloat32'
        precision = 'float32'; bytes = 4; components = 2;
    otherwise
        error('CRSW:Input', 'Unknown raster format: %s.', format);
end
fid = fopen(filename, 'r', 'ieee-be');
assert(fid >= 0, 'CRSW:Input', 'Cannot open %s.', filename);
cleanup = onCleanup(@() fclose(fid));
rows = bounds(2) - bounds(1) + 1;
cols = bounds(4) - bounds(3) + 1;
data = complex(zeros(rows, cols, 'single'));
for r = 1:rows
    offset = ((bounds(1) + r - 2) * width + bounds(3) - 1) * components * bytes;
    assert(fseek(fid, offset, 'bof') == 0, 'CRSW:Input', 'Cannot seek in %s.', filename);
    [values, count] = fread(fid, components * cols, [precision, '=>single']);
    assert(count == components * cols, 'CRSW:Input', 'Incomplete raster: %s.', filename);
    data(r,:) = complex(values(1:2:end), values(2:2:end)).';
end
end
