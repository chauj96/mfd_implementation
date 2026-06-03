function writeExtrudedMeshVTP(filename, V3, cell_struct, face_struct, cellDataStruct)
% Write extruded 3D mesh after cell classification
% The fifth argument 'cellDataStruct' MUST be a struct mapping field names to
% numeric vectors of length nCells. Each struct field will be exported as a
% separate <DataArray> inside the <CellData> block with the field name used
% as the DataArray Name.
%
% Example (recommended usage):
%  writeExtrudedMeshVTP(fname, V3, cell_struct, face_struct, ...
%      struct('saturation', s_cell_data, 'pressure', p_cell_data))
%
% Notes:
%  - All vectors must have length equal to the number of cells (nCells).
%  - Scalars will be expanded to the full cell length.
%  - Integer-like arrays are written with VTK type Int32; otherwise Float64.
%  - This simplified interface intentionally supports only struct input so
%    users explicitly provide the set of cell fields they want to visualize.

    % Validate input count
    if nargin < 5
        error('writeExtrudedMeshVTP requires five arguments; the fifth must be a struct of cell fields.');
    end

    nCells = numel(cell_struct);
    nPts   = size(V3,1);
    VTK_POLYHEDRON = 42;

    connectivity = [];
    offsets = zeros(nCells,1);
    types = VTK_POLYHEDRON * ones(nCells,1,'uint8');

    % faces/faceoffsets: per cell
    faces = [];
    faceoffsets = zeros(nCells,1);

    off_conn = 0;
    off_face = 0;

    for c = 1:nCells
        fids = cell_struct(c).faces(:)';

        % points used by the polyhedron cell = union of face vertices
        vids = [];
        for k = 1:numel(fids)
            vids = [vids, face_struct(fids(k)).verts(:)'];
        end
        vids = unique(vids, 'stable'); % 1-based
        vids0 = vids - 1; % 0-based for VTK

        connectivity = [connectivity, vids0];
        off_conn = off_conn + numel(vids0);
        offsets(c) = off_conn;

        % faces encoding for this cell:
        % [numFaces, nV(f1), v..., nV(f2), v..., ...]
        rec = numel(fids);
        for k = 1:numel(fids)
            v = face_struct(fids(k)).verts(:)' - 1; % 0-based
            rec = [rec, numel(v), v];
        end

        faces = [faces, rec];
        off_face = off_face + numel(rec);
        faceoffsets(c) = off_face;
    end

    % Prepare cell data entries: accept only a struct with named fields
    if ~isstruct(cellDataStruct)
        error(['writeExtrudedMeshVTP: this version only accepts a struct for cell data.\n', ...
               'Pass data as struct(''saturation'', s_cell_data, ''pressure'', p_cell_data).']);
    end

    fn = fieldnames(cellDataStruct);
    dataNames = cell(1, numel(fn));
    dataValues = cell(1, numel(fn));
    for k = 1:numel(fn)
        dataNames{k} = fn{k};
        dataValues{k} = cellDataStruct.(fn{k});
    end

    % Validate and normalize vectors
    for k = 1:numel(dataValues)
        v = dataValues{k};
        if isscalar(v)
            dataValues{k} = repmat(v, nCells, 1);
        else
            v = v(:);
            if numel(v) ~= nCells
                error('Cell data ''%s'' length mismatch: expected %d, got %d', dataNames{k}, nCells, numel(v));
            end
            dataValues{k} = v;
        end
    end

    fid = fopen(filename,'w');

    fprintf(fid,'<?xml version="1.0"?>\n');
    fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid,'<UnstructuredGrid>\n');
    fprintf(fid,'<Piece NumberOfPoints="%d" NumberOfCells="%d">\n', nPts, nCells);

    % Points
    fprintf(fid,'<Points>\n');
    fprintf(fid,'<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid,'%.15g %.15g %.15g\n', V3');
    fprintf(fid,'</DataArray>\n');
    fprintf(fid,'</Points>\n');

    % Cells
    fprintf(fid,'<Cells>\n');

    fprintf(fid,'<DataArray type="Int32" Name="connectivity" format="ascii">\n');
    fprintf(fid,'%d ', connectivity);
    fprintf(fid,'\n</DataArray>\n');

    fprintf(fid,'<DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fid,'%d ', offsets);
    fprintf(fid,'\n</DataArray>\n');

    fprintf(fid,'<DataArray type="UInt8" Name="types" format="ascii">\n');
    fprintf(fid,'%d ', types);
    fprintf(fid,'\n</DataArray>\n');

    % Polyhedron-specific arrays
    fprintf(fid,'<DataArray type="Int32" Name="faces" format="ascii">\n');
    fprintf(fid,'%d ', faces);
    fprintf(fid,'\n</DataArray>\n');

    fprintf(fid,'<DataArray type="Int32" Name="faceoffsets" format="ascii">\n');
    fprintf(fid,'%d ', faceoffsets);
    fprintf(fid,'\n</DataArray>\n');

    fprintf(fid,'</Cells>\n');

    % CellData: write all provided fields
    fprintf(fid,'<CellData>\n');
    for k = 1:numel(dataNames)
        name = dataNames{k};
        arr  = dataValues{k}(:);

        % Determine VTK data type: use Int32 for integer-like data, Float64 otherwise
        if all(arr == floor(arr))
            vtkType = 'Int32';
            fmt = '%d ';
            % Ensure integer type for printing
            arr_print = int32(arr);
        else
            vtkType = 'Float64';
            fmt = '%.15g ';
            arr_print = arr;
        end

        fprintf(fid,'<DataArray type="%s" Name="%s" format="ascii">\n', vtkType, name);
        fprintf(fid, fmt, arr_print);
        fprintf(fid,'\n</DataArray>\n');
    end
    fprintf(fid,'</CellData>\n');

    fprintf(fid,'</Piece>\n');
    fprintf(fid,'</UnstructuredGrid>\n');
    fprintf(fid,'</VTKFile>\n');

    fclose(fid);
    fprintf('Wrote %s\n', filename);
end
