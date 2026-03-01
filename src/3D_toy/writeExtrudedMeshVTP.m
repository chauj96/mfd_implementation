function writeExtrudedMeshVTP(filename, V3, cell_struct, face_struct, cellData, cellDataName)
% Write extruded 3D mesh after cell classification
    
    if nargin < 6, cellDataName = 'cellMarking'; end
    
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
    
    % CellData
    fprintf(fid,'<CellData Scalars="%s">\n', cellDataName);
    fprintf(fid,'<DataArray type="Int32" Name="%s" format="ascii">\n', cellDataName);
    fprintf(fid,'%d\n', cellData);
    fprintf(fid,'</DataArray>\n');
    fprintf(fid,'</CellData>\n');
    
    fprintf(fid,'</Piece>\n');
    fprintf(fid,'</UnstructuredGrid>\n');
    fprintf(fid,'</VTKFile>\n');
    
    fclose(fid);
    fprintf('Wrote %s\n', filename);
end
