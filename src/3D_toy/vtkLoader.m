function [cell_struct, face_struct, V3, cells3D, attribute] = vtkLoader(filename)
% vtkLoader  Load a binary VTU (UnstructuredGrid) file and build the same
%            cell_struct / face_struct arrays that MRSTGridConvert produces.
%
%   Supports:
%     - VTK type 10  (tetrahedra,   4 nodes, 4 triangular faces)
%     - VTK type 12  (hexahedra,    8 nodes, 6 quad faces)
%     - VTK type 13  (wedges,       6 nodes, 2 tri + 3 quad faces)
%     - VTK type 14  (pyramids,     5 nodes, 1 quad + 4 tri faces)
%     - VTK type 42  (polyhedra,    arbitrary faces via face_connectivity etc.)
%     - mixed meshes combining any of the above
%     - binary format with vtkZLibDataCompressor
%     - header_type="UInt32" (4-byte block headers) and "UInt64" (8-byte)

    %% ------------------------------------------------------------------ %%
    %% 1.  Read XML header & locate DataArray blocks
    %% ------------------------------------------------------------------ %%
    fid = fopen(filename, 'rb');
    if fid < 0, error('vtkLoader: cannot open file "%s"', filename); end
    raw = fread(fid, Inf, '*uint8');
    fclose(fid);
    xml_str = char(raw(:)');

    % Parse header_type (UInt32 = 4-byte block-size words, UInt64 = 8-byte)
    hdr_type_tok = regexpi(xml_str, 'header_type\s*=\s*[''"]([^''"]*)[''""]', 'tokens', 'once');
    if ~isempty(hdr_type_tok) && strcmpi(hdr_type_tok{1}, 'UInt64')
        hdr_word_bytes = 8;
    else
        hdr_word_bytes = 4;  % default / UInt32
    end

    piece_tok = regexpi(xml_str, '<Piece\s+([^>]*)>', 'tokens', 'once');
    if isempty(piece_tok), error('vtkLoader: no <Piece> tag found'); end
    piece_str = piece_tok{1};
    nPoints = str2double(xmlAttr(piece_str, 'NumberOfPoints'));
    nCells  = str2double(xmlAttr(piece_str, 'NumberOfCells'));
    fprintf('vtkLoader: %d points, %d cells\n', nPoints, nCells);

    %% ------------------------------------------------------------------ %%
    %% 2.  Decode all DataArray blocks
    %% ------------------------------------------------------------------ %%
    [da_tags, da_starts, da_ends] = findDataArrays(xml_str);

    pts_data  = getDataByName(da_tags, da_starts, da_ends, xml_str, 'Points',              hdr_word_bytes);
    conn_data = getDataByName(da_tags, da_starts, da_ends, xml_str, 'connectivity',        hdr_word_bytes);
    off_data  = getDataByName(da_tags, da_starts, da_ends, xml_str, 'offsets',             hdr_word_bytes);
    typ_data  = getDataByName(da_tags, da_starts, da_ends, xml_str, 'types',               hdr_word_bytes);

    % Polyhedral-specific arrays (may be absent for non-poly meshes)
    fconn_data  = getDataByName(da_tags, da_starts, da_ends, xml_str, 'face_connectivity', hdr_word_bytes);
    foff_data   = getDataByName(da_tags, da_starts, da_ends, xml_str, 'face_offsets',      hdr_word_bytes);
    p2f_data    = getDataByName(da_tags, da_starts, da_ends, xml_str, 'polyhedron_to_faces', hdr_word_bytes);
    poff_data   = getDataByName(da_tags, da_starts, da_ends, xml_str, 'polyhedron_offsets', hdr_word_bytes);

    % Cell-data attribute (try common names)
    attr_data = getDataByName(da_tags, da_starts, da_ends, xml_str, 'attribute',           hdr_word_bytes);
    if isempty(attr_data)
        attr_data = getDataByName(da_tags, da_starts, da_ends, xml_str, 'CellEntityIds',   hdr_word_bytes);
    end

    %% ------------------------------------------------------------------ %%
    %% 3.  Cast raw bytes to typed arrays
    %%     Use the 'type' attribute from the DataArray XML tag to pick
    %%     the correct MATLAB typecast target.
    %% ------------------------------------------------------------------ %%

    % Points: Float32 or Float64
    pts_type = getTypeByName(da_tags, 'Points');
    if strcmpi(pts_type, 'Float64')
        V3 = double(reshape(typecast(pts_data, 'double'), 3, [])');
    else
        V3 = double(reshape(typecast(pts_data, 'single'), 3, [])');
    end

    % offsets: Int32 or Int64
    off_type = getTypeByName(da_tags, 'offsets');
    if strcmpi(off_type, 'Int32')
        offs = double(typecast(off_data, 'int32'));
    else
        offs = double(typecast(off_data, 'int64'));
    end

    typs = typecast(typ_data, 'uint8');                              % nCells x 1

    if ~isempty(attr_data)
        % Support both Int32 and Int64 attribute arrays
        if numel(attr_data) == nCells * 4
            attribute = double(typecast(attr_data, 'int32'));
        else
            attribute = double(typecast(attr_data, 'int64'));
        end
    else
        attribute = [];
    end

    % Standard connectivity: Int32 or Int64
    conn_type = getTypeByName(da_tags, 'connectivity');
    if strcmpi(conn_type, 'Int32')
        conn = double(reshape(typecast(conn_data, 'int32'), 1, [])) + 1; % 1-based
    else
        conn = double(reshape(typecast(conn_data, 'int64'), 1, [])) + 1; % 1-based
    end

    % Check for type 42
    hasPolyhedra = any(typs == 42);
    supportedTypes = [10, 12, 13, 14, 42];
    uniqueTypes = unique(typs);
    unsupported = setdiff(uniqueTypes, supportedTypes);
    if ~isempty(unsupported)
        warning('vtkLoader: unsupported cell types: %s', num2str(unsupported'));
    end

    %% ------------------------------------------------------------------ %%
    %% 4a.  Decode polyhedral face tables (type 42 only)
    %%
    %%  face_connectivity : flat node list for ALL faces (Int64, 0-based)
    %%  face_offsets      : face_offsets(i) = end index in face_connectivity
    %%                      for face i  (length = total number of named faces)
    %%  polyhedron_to_faces : flat list of face indices (0-based) for all
    %%                      polyhedral cells concatenated
    %%  polyhedron_offsets  : polyhedron_offsets(ci) = end index into
    %%                      polyhedron_to_faces for polyhedral cell ci
    %%                      (length = number of type-42 cells)
    %% ------------------------------------------------------------------ %%
    poly_cell_faces = {};   % poly_cell_faces{k} = face indices (1-based global) for poly cell k
    poly_face_verts = {};   % poly_face_verts{f} = node list (1-based) for global face f

    if hasPolyhedra
        if isempty(fconn_data) || isempty(foff_data) || isempty(p2f_data) || isempty(poff_data)
            error('vtkLoader: type-42 cells detected but face_connectivity/face_offsets/polyhedron_to_faces/polyhedron_offsets arrays are missing');
        end

        fconn_type = getTypeByName(da_tags, 'face_connectivity');
        if strcmpi(fconn_type, 'Int32')
            fc_nodes = double(reshape(typecast(fconn_data, 'int32'), 1, [])) + 1;
        else
            fc_nodes = double(reshape(typecast(fconn_data, 'int64'), 1, [])) + 1;
        end

        foff_type = getTypeByName(da_tags, 'face_offsets');
        if strcmpi(foff_type, 'Int32')
            fc_offs = double(typecast(foff_data, 'int32'));
        else
            fc_offs = double(typecast(foff_data, 'int64'));
        end

        p2f_type = getTypeByName(da_tags, 'polyhedron_to_faces');
        if strcmpi(p2f_type, 'Int32')
            p2f_flat = double(typecast(p2f_data, 'int32')) + 1;
        else
            p2f_flat = double(typecast(p2f_data, 'int64')) + 1;
        end

        poff_type = getTypeByName(da_tags, 'polyhedron_offsets');
        if strcmpi(poff_type, 'Int32')
            p_offs = double(typecast(poff_data, 'int32'));
        else
            p_offs = double(typecast(poff_data, 'int64'));
        end

        nPolyFaces = numel(fc_offs);
        poly_face_verts = cell(nPolyFaces, 1);
        prev = 0;
        for pf = 1:nPolyFaces
            poly_face_verts{pf} = fc_nodes(prev+1 : fc_offs(pf));
            prev = fc_offs(pf);
        end

        % Map polyhedral cells (in order they appear in type array) to face lists
        nPolyCells = numel(p_offs);
        poly_cell_faces = cell(nPolyCells, 1);
        prev = 0;
        for pk = 1:nPolyCells
            poly_cell_faces{pk} = p2f_flat(prev+1 : p_offs(pk));
            prev = p_offs(pk);
        end
    end

    %% ------------------------------------------------------------------ %%
    %% 4b.  Build per-cell connectivity for standard (non-42) cells
    %% ------------------------------------------------------------------ %%
    % Local face tables (1-based local node indices per cell type)
    localFaceDef = cell(15,1);
    localFaceDef{10} = {[1 2 4]; [2 3 4]; [3 1 4]; [1 3 2]};
    localFaceDef{12} = {[1 4 3 2]; [5 6 7 8]; [1 2 6 5]; [2 3 7 6]; [3 4 8 7]; [4 1 5 8]};
    localFaceDef{13} = {[1 2 3]; [4 6 5]; [1 4 5 2]; [2 5 6 3]; [3 6 4 1]};
    localFaceDef{14} = {[1 4 3 2]; [1 2 5]; [2 3 5]; [3 4 5]; [4 1 5]};

    cell_nodes = cell(nCells, 1);
    prev = 0;
    polyIdx = 0;   % counter into poly_cell_faces
    for ci = 1:nCells
        tp = typs(ci);
        if tp == 42
            polyIdx = polyIdx + 1;
            % Collect unique vertices from all faces of this poly cell
            fids_local = poly_cell_faces{polyIdx};
            all_verts = [];
            for k = 1:numel(fids_local)
                all_verts = [all_verts, poly_face_verts{fids_local(k)}]; %#ok<AGROW>
            end
            cell_nodes{ci} = unique(all_verts(:))';
            prev = offs(ci);  % advance past connectivity entry (unused for poly)
        else
            cell_nodes{ci} = conn(prev+1 : offs(ci));
            prev = offs(ci);
        end
    end

    %% ------------------------------------------------------------------ %%
    %% 5.  Build half-face list from all cells
    %%
    %%  For type-42 cells: half-faces come directly from poly_face_verts.
    %%  For other types:   half-faces come from localFaceDef.
    %%
    %%  Key for face matching: sorted node ids (padded to max face size).
    %% ------------------------------------------------------------------ %%

    % Determine max face size for key width
    maxFaceNodes = 4;  % quads for standard types
    if hasPolyhedra
        for pf = 1:numel(poly_face_verts)
            maxFaceNodes = max(maxFaceNodes, numel(poly_face_verts{pf}));
        end
    end

    % Count total half-faces
    nHF = 0;
    polyIdx = 0;
    for ci = 1:nCells
        tp = typs(ci);
        if tp == 42
            polyIdx = polyIdx + 1;
            nHF = nHF + numel(poly_cell_faces{polyIdx});
        else
            nHF = nHF + numel(localFaceDef{tp});
        end
    end

    HF_key   = zeros(nHF, maxFaceNodes, 'int32');
    HF_verts = cell(nHF, 1);
    HF_cell  = zeros(nHF, 1, 'int32');

    idx = 0;
    polyIdx = 0;
    for ci = 1:nCells
        tp = typs(ci);
        if tp == 42
            polyIdx = polyIdx + 1;
            fids_local = poly_cell_faces{polyIdx};
            for k = 1:numel(fids_local)
                idx = idx + 1;
                vv = poly_face_verts{fids_local(k)};
                sv = sort(double(vv));
                HF_key(idx, 1:numel(sv)) = int32(sv);
                HF_verts{idx} = vv;
                HF_cell(idx)  = ci;
            end
        else
            cnodes = cell_nodes{ci};
            lfd    = localFaceDef{tp};
            for lf = 1:numel(lfd)
                idx = idx + 1;
                vv  = cnodes(lfd{lf});
                sv  = sort(double(vv));
                HF_key(idx, 1:numel(sv)) = int32(sv);
                HF_verts{idx} = vv;
                HF_cell(idx)  = ci;
            end
        end
    end

    %% ------------------------------------------------------------------ %%
    %% 5b.  Match half-faces to unique faces
    %% ------------------------------------------------------------------ %%
    [HF_key_s, ia] = sortrows(HF_key);
    HF_verts_s = HF_verts(ia);
    HF_cell_s  = HF_cell(ia);

    nFaces       = 0;
    face_verts_c = cell(nHF, 1);
    face_cells   = zeros(nHF, 2, 'int32');
    i = 1;
    while i <= nHF
        j = i;
        while j < nHF && all(HF_key_s(j+1,:) == HF_key_s(i,:))
            j = j + 1;
        end
        nFaces = nFaces + 1;
        face_verts_c{nFaces} = HF_verts_s{i};
        if j == i
            face_cells(nFaces,1) = HF_cell_s(i);
            face_cells(nFaces,2) = 0;
        else
            face_cells(nFaces,1) = HF_cell_s(i);
            face_cells(nFaces,2) = HF_cell_s(i+1);
        end
        i = j + 1;
    end
    face_verts_c = face_verts_c(1:nFaces);
    face_cells   = face_cells(1:nFaces,:);
    fprintf('vtkLoader: %d unique faces\n', nFaces);

    %% ------------------------------------------------------------------ %%
    %% 6.  Cell -> face map
    %% ------------------------------------------------------------------ %%
    cell_faces_list = cell(nCells, 1);
    for f = 1:nFaces
        c1 = face_cells(f,1);  if c1 > 0, cell_faces_list{c1} = [cell_faces_list{c1}, f]; end
        c2 = face_cells(f,2);  if c2 > 0, cell_faces_list{c2} = [cell_faces_list{c2}, f]; end
    end

    %% ------------------------------------------------------------------ %%
    %% 7.  Face geometry
    %%     Fan-triangulate each polygon from its first vertex.
    %%     Centroid = area-weighted average of triangle centroids.
    %%     Normal   = sum of triangle cross-products (area-weighted direction).
    %% ------------------------------------------------------------------ %%
    face_centroid = zeros(nFaces, 3);
    face_normal   = zeros(nFaces, 3);
    face_area     = zeros(nFaces, 1);

    for f = 1:nFaces
        vids = double(face_verts_c{f});
        pts  = V3(vids, :);
        nv   = size(pts, 1);
        if nv < 3
            face_centroid(f,:) = mean(pts, 1);
            continue;
        end
        cr_sum  = [0 0 0];
        xc_sum  = [0 0 0];
        area_tot = 0;
        p0 = pts(1,:);
        for kv = 2:nv-1
            cr = cross(pts(kv,:) - p0, pts(kv+1,:) - p0);
            tri_area2 = norm(cr);
            cr_sum   = cr_sum  + cr;
            xc_sum   = xc_sum  + tri_area2 * (p0 + pts(kv,:) + pts(kv+1,:));
            area_tot = area_tot + tri_area2;
        end
        total_area = area_tot / 2;
        face_area(f) = total_area;
        if area_tot > 0
            face_normal(f,:)   = cr_sum / norm(cr_sum);   % unit normal
            face_centroid(f,:) = xc_sum / (3 * area_tot); % area-weighted centroid
        else
            face_centroid(f,:) = mean(pts, 1);
        end
    end

    %% ------------------------------------------------------------------ %%
    %% 8.  Cell geometry (centroid & volume)
    %%     Works for any cell type by decomposing faces into tetrahedra
    %%     from an approximate centroid (vertex mean of the cell).
    %%
    %%     V = sum_f sum_tri  dot(a, cross(b,c)) / 6
    %%     xC = sum_f sum_tri sv * (xc_approx + p0 + pk + pk+1) / 4  / V
    %% ------------------------------------------------------------------ %%
    cell_centroid = zeros(nCells, 3);
    cell_volume   = zeros(nCells, 1);

    polyIdx = 0;
    for ci = 1:nCells
        tp = typs(ci);

        % Build list of face vertex lists for this cell
        if tp == 42
            polyIdx = polyIdx + 1;
            fids_local = poly_cell_faces{polyIdx};
            cell_face_verts = cell(numel(fids_local), 1);
            for k = 1:numel(fids_local)
                cell_face_verts{k} = double(poly_face_verts{fids_local(k)});
            end
        else
            cnodes = double(cell_nodes{ci});
            lfd    = localFaceDef{tp};
            cell_face_verts = cell(numel(lfd), 1);
            for lf = 1:numel(lfd)
                cell_face_verts{lf} = cnodes(lfd{lf});
            end
        end

        % Approximate centroid: average of unique vertex positions
        all_v = unique(cat(2, cell_face_verts{:}));
        xc_approx = mean(V3(all_v, :), 1);

        vol  = 0;
        xc_w = [0 0 0];
        for lf = 1:numel(cell_face_verts)
            fvids = cell_face_verts{lf};
            fpts  = V3(fvids, :);
            nfv   = size(fpts, 1);
            for kv = 2:nfv-1
                a  = fpts(1,:)    - xc_approx;
                b  = fpts(kv,:)   - xc_approx;
                c  = fpts(kv+1,:) - xc_approx;
                sv = dot(a, cross(b, c)) / 6;
                vol  = vol  + sv;
                % tet centroid = (xc_approx + p0 + pk + pk+1) / 4
                xc_w = xc_w + sv * (xc_approx + fpts(1,:) + fpts(kv,:) + fpts(kv+1,:)) / 4;
            end
        end
        cell_volume(ci) = abs(vol);
        if abs(vol) > 0
            % xc_w / vol  (use signed vol so direction is consistent)
            cell_centroid(ci,:) = xc_w / vol;
        else
            cell_centroid(ci,:) = xc_approx;
        end
    end

    %% ------------------------------------------------------------------ %%
    %% 9.  Populate face_struct
    %% ------------------------------------------------------------------ %%
    face_struct = struct('cells',  cell(1,nFaces), ...
                         'verts',  cell(1,nFaces), ...
                         'center', cell(1,nFaces), ...
                         'normal', cell(1,nFaces), ...
                         'area',   cell(1,nFaces));
    for f = 1:nFaces
        neigh = face_cells(f,:);
        neigh = neigh(neigh > 0);
        face_struct(f).cells  = neigh(:)';
        face_struct(f).verts  = int32(face_verts_c{f}(:));
        face_struct(f).center = face_centroid(f,:)';
        face_struct(f).normal = face_normal(f,:)';
        face_struct(f).area   = face_area(f);
    end

    %% ------------------------------------------------------------------ %%
    %% 10. Populate cell_struct
    %% ------------------------------------------------------------------ %%
    cell_struct = struct('center',            cell(1,nCells), ...
                         'volume',            cell(1,nCells), ...
                         'faces',             cell(1,nCells), ...
                         'faces_orientation', cell(1,nCells), ...
                         'face_normals',      cell(1,nCells));
    cells3D = cell(nCells, 1);

    for ci = 1:nCells
        fids = cell_faces_list{ci};
        nf   = numel(fids);
        cell_struct(ci).center = cell_centroid(ci,:)';
        cell_struct(ci).volume = cell_volume(ci);
        cell_struct(ci).faces  = fids;

        face_signs   = zeros(1, nf);
        face_normals = zeros(nf, 3);
        xc = cell_centroid(ci,:)';
        for k = 1:nf
            f      = fids(k);
            nf_vec = face_struct(f).normal(:);
            xf     = face_struct(f).center(:);
            d      = dot(nf_vec, xf - xc);
            if d < 0
                face_signs(k)     = -1;
                face_normals(k,:) = -nf_vec';
            else
                face_signs(k)     =  1;
                face_normals(k,:) =  nf_vec';
            end
        end
        cell_struct(ci).faces_orientation = face_signs;
        cell_struct(ci).face_normals      = face_normals;
        cells3D{ci} = unique(double(cell_nodes{ci}(:)), 'stable')';
    end

    fprintf('vtkLoader: done.\n');
    fprintf('  %d vertices, %d cells, %d faces\n', size(V3,1), nCells, nFaces);
end

%% ======================================================================
%% LOCAL HELPERS
%% ======================================================================

function val = xmlAttr(tag, attr)
    tok = regexpi(tag, [attr '\s*=\s*[''"]([^''"]*)[''"]'], 'tokens', 'once');
    if isempty(tok), val = ''; else, val = tok{1}; end
end

function [da_tags, da_starts, da_ends] = findDataArrays(xml_str)
    starts = regexp(xml_str, '<DataArray\s', 'start');
    da_tags   = cell(numel(starts),1);
    da_starts = zeros(numel(starts),1);
    da_ends   = zeros(numel(starts),1);
    for k = 1:numel(starts)
        tag_end = find(xml_str(starts(k):end) == '>', 1) + starts(k) - 1;
        da_tags{k}   = xml_str(starts(k):tag_end);
        da_starts(k) = tag_end + 1;
        close_pos    = strfind(xml_str(tag_end:end), '</DataArray>');
        da_ends(k)   = tag_end + close_pos(1) - 2;
    end
end

function raw_bytes = getDataByName(da_tags, da_starts, da_ends, xml_str, name, hdr_word_bytes)
    raw_bytes = [];
    for k = 1:numel(da_tags)
        tok = regexpi(da_tags{k}, 'Name\s*=\s*[''"]([^''"]*)[''""]', 'tokens', 'once');
        if isempty(tok), continue; end
        if ~strcmpi(tok{1}, name), continue; end
        b64 = strtrim(xml_str(da_starts(k):da_ends(k)));
        lt_pos = find(b64 == '<', 1);
        if ~isempty(lt_pos), b64 = strtrim(b64(1:lt_pos-1)); end
        if isempty(b64), continue; end
        raw_bytes = decodeVTKBinary(b64, hdr_word_bytes);
        return;
    end
end

function tp = getTypeByName(da_tags, name)
% Return the 'type' attribute (e.g. 'Float32', 'Int64') for a DataArray by Name.
% Returns '' if not found; caller should default to Int64 / Float32.
    tp = '';
    for k = 1:numel(da_tags)
        tok = regexpi(da_tags{k}, 'Name\s*=\s*[''"]([^''"]*)[''""]', 'tokens', 'once');
        if isempty(tok), continue; end
        if ~strcmpi(tok{1}, name), continue; end
        type_tok = regexpi(da_tags{k}, '\btype\s*=\s*[''"]([^''"]*)[''""]', 'tokens', 'once');
        if ~isempty(type_tok), tp = type_tok{1}; end
        return;
    end
end

function data = decodeVTKBinary(b64str, hdr_word_bytes)
% Decode a VTK binary DataArray encoded with vtkZLibDataCompressor.
%
% VTK compressed binary format (little-endian, each word = hdr_word_bytes):
%   word 0          : nblocks
%   word 1          : uncompressed block size
%   word 2          : uncompressed size of last block
%   words 3..2+n    : compressed byte-size of each block (n = nblocks)
%   [compressed block 0] ... [compressed block n-1]
%
% The header and each compressed block are each independently base64-
% encoded and concatenated.  '=' padding may appear between chunks.

    if nargin < 2 || isempty(hdr_word_bytes)
        hdr_word_bytes = 4;
    end

    % Strip whitespace/newlines
    b64str = regexprep(b64str, '\s+', '');
    if isempty(b64str)
        error('vtkLoader:decodeVTKBinary', 'Empty base64 string');
    end

    % Helper to read uint word from raw bytes
    function v = readWord(raw, byteIdx)
        if hdr_word_bytes == 4
            v = double(typecast(raw(byteIdx:byteIdx+3), 'uint32'));
        else
            v = double(typecast(raw(byteIdx:byteIdx+7), 'uint64'));
        end
    end

    % ---- Split concatenated base64 string into individual chunks ----
    % Each chunk is independently padded to a multiple of 4 characters.
    % We split on boundaries where '=' padding ends and a new chunk begins.
    % Strategy: scan for every group of 4 characters.  If a group contains
    % '=', it is the last group of a chunk.
    chunks = splitB64Chunks(b64str);

    % ---- Decode header chunk (always the first chunk) ----
    hdr_raw = b64decode_one(chunks{1});

    % Probe nblocks from first word
    if numel(hdr_raw) < hdr_word_bytes
        error('vtkLoader:decodeVTKBinary', 'Header too short');
    end
    nblocks = readWord(hdr_raw, 1);
    if ~(nblocks >= 1 && nblocks <= 1e7)
        error('vtkLoader:decodeVTKBinary', 'Implausible nblocks=%g', nblocks);
    end
    header_bytes = (3 + nblocks) * hdr_word_bytes;

    % The header may span more than one chunk if the first chunk
    % didn't contain enough bytes.  Usually it's exactly one chunk.
    % If we got fewer bytes than needed, try monolithic decode.
    if numel(hdr_raw) < header_bytes
        % Fall back to decoding all chunks joined without padding
        all_raw = b64decode_one(strjoin(chunks, ''));
        hdr_raw = all_raw(1:header_bytes);
        nblocks = readWord(hdr_raw, 1);
        header_bytes = (3 + nblocks) * hdr_word_bytes;
    end
    hdr_raw = hdr_raw(1:header_bytes);

    % Read compressed sizes from header
    csizes = zeros(1, nblocks, 'double');
    for b = 1:nblocks
        off = 3*hdr_word_bytes + (b-1)*hdr_word_bytes + 1;
        csizes(b) = readWord(hdr_raw, off);
    end

    % ---- Decode compressed blocks ----
    % Case 1: chunked layout — we have exactly (1 + nblocks) chunks
    nChunks = numel(chunks);
    if nChunks == 1 + nblocks
        blks = cell(nblocks, 1);
        for b = 1:nblocks
            rb = b64decode_one(chunks{1 + b});
            cs = csizes(b);
            if numel(rb) < cs
                error('vtkLoader:decodeVTKBinary', ...
                    'Block %d: decoded %d bytes but expected %d', b, numel(rb), cs);
            end
            blks{b} = zlib_decompress(rb(1:cs));
        end
        data = vertcat(blks{:});
    elseif nChunks == 1
        % Case 2: monolithic layout — everything in one base64 blob
        if ~exist('all_raw','var')
            all_raw = hdr_raw;  % Already decoded above, but may be partial
            % Re-decode the whole thing
            all_raw = b64decode_one(chunks{1});
        end
        pos = header_bytes + 1;
        blks = cell(nblocks, 1);
        for b = 1:nblocks
            cs = csizes(b);
            if pos + cs - 1 > numel(all_raw)
                error('vtkLoader:decodeVTKBinary', ...
                    'Block %d overruns buffer (pos=%d, cs=%d, len=%d)', ...
                    b, pos, cs, numel(all_raw));
            end
            blks{b} = zlib_decompress(all_raw(pos:pos+cs-1));
            pos = pos + cs;
        end
        data = vertcat(blks{:});
    else
        % Case 3: chunk count doesn't match expected — maybe some blocks
        % share a chunk, or VTK used a different grouping.
        % Concatenate all chunks after the header and decode as one blob.
        rest_str = strjoin(chunks(2:end), '');
        rest_raw = b64decode_one(rest_str);
        pos = 1;
        blks = cell(nblocks, 1);
        for b = 1:nblocks
            cs = csizes(b);
            if pos + cs - 1 > numel(rest_raw)
                error('vtkLoader:decodeVTKBinary', ...
                    'Block %d overruns rest buffer (pos=%d, cs=%d, len=%d)', ...
                    b, pos, cs, numel(rest_raw));
            end
            blks{b} = zlib_decompress(rest_raw(pos:pos+cs-1));
            pos = pos + cs;
        end
        data = vertcat(blks{:});
    end

    data = data(:);
end


function chunks = splitB64Chunks(s)
% Split a concatenated base64 string into individual chunks.
% Each independently-encoded chunk is padded to a multiple of 4 chars.
% We detect boundaries by finding '=' padding followed by a non-'=' char.
    n = numel(s);
    chunks = {};
    start = 1;
    i = 1;
    while i <= n
        if s(i) == '='
            % Scan past all '=' chars (there can be 1 or 2)
            j = i;
            while j <= n && s(j) == '='
                j = j + 1;
            end
            % End of a chunk at position j-1
            chunks{end+1} = s(start:j-1); %#ok<AGROW>
            start = j;
            i = j;
        else
            i = i + 1;
        end
    end
    % Remaining characters (last chunk may have no padding)
    if start <= n
        chunks{end+1} = s(start:n);
    end
    if isempty(chunks)
        chunks = {s};
    end
end


function out = b64decode_one(str)
% Decode a single base64-encoded chunk to uint8 bytes.
% Strips any internal '=' that are NOT at the very end (for joined strings),
% and ensures proper end-padding.
    str = regexprep(str, '\s+', '');
    if isempty(str)
        out = uint8([]);
        return;
    end

    % Remove all '=' first, then re-pad at the end
    str_clean = strrep(str, '=', '');
    if isempty(str_clean)
        out = uint8([]);
        return;
    end
    rem4 = mod(numel(str_clean), 4);
    if rem4 == 2
        str_clean = [str_clean '=='];
    elseif rem4 == 3
        str_clean = [str_clean '='];
    elseif rem4 == 1
        % Drop last char (invalid)
        str_clean = str_clean(1:end-1);
        rem4 = mod(numel(str_clean), 4);
        if rem4 == 2, str_clean = [str_clean '==']; end
        if rem4 == 3, str_clean = [str_clean '=']; end
    end

    decoder = java.util.Base64.getDecoder();
    jbytes  = decoder.decode(str_clean);
    out     = typecast(jbytes, 'uint8');
    out     = out(:);
end

function out = zlib_decompress(data)
    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier();
    try
        a = java.io.ByteArrayInputStream(data);
        b = java.util.zip.InflaterInputStream(a);
        c = java.io.ByteArrayOutputStream();
        isc.copyStream(b, c);
        out = typecast(c.toByteArray(), 'uint8');
        out = out(:);
    catch
        try
            inf = java.util.zip.Inflater(true);
            inf.setInput(data);
            buf = java.io.ByteArrayOutputStream();
            tmp = zeros(1, 65536, 'int8');
            while ~inf.finished()
                n = inf.inflate(tmp);
                if n > 0, buf.write(tmp, 0, n); end
                if n == 0 && ~inf.finished(), break; end
            end
            inf.end();
            out = typecast(buf.toByteArray(), 'uint8');
            out = out(:);
        catch ME2
            error('vtkLoader:zlib_decompress', ...
                'Could not decompress: %s', ME2.message);
        end
    end
end
