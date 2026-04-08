function [cell_struct, face_struct, V3, cells3D, attribute] = vtkLoader(filename)
% vtkLoader  Load a binary VTU (UnstructuredGrid) file and build the same
%            cell_struct / face_struct arrays that MRSTGridConvert produces.
%
%   [cell_struct, face_struct, V3, cells3D, attribute] = vtkLoader(filename)
%
%   Supports:
%     - VTK type 10  (tetrahedra,   4 nodes, 4 triangular faces)
%     - VTK type 12  (hexahedra,    8 nodes, 6 quad faces)
%     - VTK type 13  (wedges,       6 nodes, 2 tri + 3 quad faces)
%     - VTK type 14  (pyramids,     5 nodes, 1 quad + 4 tri faces)
%     - mixed meshes combining any of the above
%     - binary format with vtkZLibDataCompressor
%     - header_type="UInt32"  (4-byte block headers)
%     - Float32 points, Int64 connectivity/offsets, UInt8 types, Int64 attribute
%
%   Outputs match MRSTGridConvert exactly:
%     V3          - nPts x 3  node coordinates (double)
%     cell_struct - 1 x nCells struct with fields:
%                     .center            (3x1)
%                     .volume            (scalar)
%                     .faces             (1 x nf  face indices, 1-based)
%                     .faces_orientation (1 x nf  ±1)
%                     .face_normals      (nf x 3)
%     face_struct - 1 x nFaces struct with fields:
%                     .cells  (1 x 1 or 1 x 2  cell indices, 1-based)
%                     .verts  (nv x 1  vertex indices, 1-based)
%                     .center (3x1)
%                     .normal (3x1  outward unit normal, sign wrt cell 1)
%                     .area   (scalar)
%     cells3D     - nCells x 1 cell array, each entry = unique vertex ids
%     attribute   - nCells x 1 int64 cell-data array (or [] if absent)

    %% ------------------------------------------------------------------ %%
    %% 1.  Read XML header & locate DataArray blocks
    %% ------------------------------------------------------------------ %%
    fid = fopen(filename, 'rb');
    if fid < 0
        error('vtkLoader: cannot open file "%s"', filename);
    end
    raw = fread(fid, Inf, '*uint8');
    fclose(fid);

    xml_str = char(raw(:)');

    % Parse Piece attributes
    piece_tok = regexpi(xml_str, '<Piece\s+([^>]*)>', 'tokens', 'once');
    if isempty(piece_tok)
        error('vtkLoader: no <Piece> tag found');
    end
    piece_str = piece_tok{1};
    nPoints = str2double(xmlAttr(piece_str, 'NumberOfPoints'));
    nCells  = str2double(xmlAttr(piece_str, 'NumberOfCells'));

    fprintf('vtkLoader: %d points, %d cells\n', nPoints, nCells);

    %% ------------------------------------------------------------------ %%
    %% 2.  Decode all DataArray blocks
    %% ------------------------------------------------------------------ %%
    % Find all DataArray tags and their base64 content
    [da_tags, da_starts, da_ends] = findDataArrays(xml_str);

    % Locate each required array by Name
    pts_data  = getDataByName(da_tags, da_starts, da_ends, xml_str, 'Points');
    conn_data = getDataByName(da_tags, da_starts, da_ends, xml_str, 'connectivity');
    off_data  = getDataByName(da_tags, da_starts, da_ends, xml_str, 'offsets');
    typ_data  = getDataByName(da_tags, da_starts, da_ends, xml_str, 'types');
    attr_data = getDataByName(da_tags, da_starts, da_ends, xml_str, 'attribute');  % may be []

    %% ------------------------------------------------------------------ %%
    %% 3.  Cast raw bytes to typed arrays
    %% ------------------------------------------------------------------ %%
    V3   = double(reshape(typecast(pts_data,  'single'), 3, [])');   % nPts x 3
    conn = double(reshape(typecast(conn_data, 'int64'),  1, [])) + 1; % 1-based
    offs = double(typecast(off_data, 'int64'));                        % nCells x 1
    typs = typecast(typ_data, 'uint8');                               % nCells x 1

    if ~isempty(attr_data)
        attribute = typecast(attr_data, 'int64');
    else
        attribute = [];
    end

    % Report cell types present
    uniqueTypes = unique(typs);
    supportedTypes = [10, 12, 13, 14];  % tet, hex, wedge, pyramid
    unsupported = setdiff(uniqueTypes, supportedTypes);
    if ~isempty(unsupported)
        warning('vtkLoader: unsupported cell types detected: %s', num2str(unsupported'));
    end

    %% ------------------------------------------------------------------ %%
    %% 4.  Build per-cell connectivity (variable nodes per cell)
    %%
    %%  VTK node ordering (0-based), converted to 1-based here:
    %%   type 10 (tet):     4 nodes
    %%   type 12 (hex):     8 nodes
    %%   type 13 (wedge):   6 nodes
    %%   type 14 (pyramid): 5 nodes
    %% ------------------------------------------------------------------ %%
    cell_nodes = cell(nCells, 1);
    prev = 0;
    for ci = 1:nCells
        cell_nodes{ci} = conn(prev+1 : offs(ci));
        prev = offs(ci);
    end

    %% ------------------------------------------------------------------ %%
    %% 5.  Build face topology from cells
    %%
    %%  Local face definitions (1-based local node indices) per VTK cell type.
    %%  Each face is a row; faces with fewer nodes are padded with 0.
    %%
    %%  type 10 (tet): 4 triangular faces
    %%    face 1: [1 2 4],  face 2: [2 3 4],  face 3: [3 1 4],  face 4: [1 3 2]
    %%
    %%  type 12 (hex): 6 quad faces (VTK standard ordering)
    %%    face 1: [1 4 3 2],  face 2: [5 6 7 8],  face 3: [1 2 6 5]
    %%    face 4: [2 3 7 6],  face 5: [3 4 8 7],  face 6: [4 1 5 8]
    %%
    %%  type 13 (wedge): 2 triangular + 3 quad faces
    %%    face 1: [1 2 3],  face 2: [4 6 5]
    %%    face 3: [1 4 5 2],  face 4: [2 5 6 3],  face 5: [3 6 4 1]
    %%
    %%  type 14 (pyramid): 1 quad base + 4 triangular faces
    %%    face 1: [1 4 3 2],  face 2: [1 2 5],  face 3: [2 3 5]
    %%    face 4: [3 4 5],    face 5: [4 1 5]
    %% ------------------------------------------------------------------ %%

    % Local face tables: cell array indexed by VTK type
    % Each entry is a cell array of face node lists (1-based local indices)
    localFaceDef = cell(15,1);
    localFaceDef{10} = {[1 2 4]; [2 3 4]; [3 1 4]; [1 3 2]};
    localFaceDef{12} = {[1 4 3 2]; [5 6 7 8]; [1 2 6 5]; [2 3 7 6]; [3 4 8 7]; [4 1 5 8]};
    localFaceDef{13} = {[1 2 3]; [4 6 5]; [1 4 5 2]; [2 5 6 3]; [3 6 4 1]};
    localFaceDef{14} = {[1 4 3 2]; [1 2 5]; [2 3 5]; [3 4 5]; [4 1 5]};

    % Count total half-faces for preallocation
    nHF = 0;
    for ci = 1:nCells
        tp = typs(ci);
        nHF = nHF + numel(localFaceDef{tp});
    end

    % Half-face storage: sorted vertex key (for matching), cell id, local face id
    % For faces with >3 nodes we sort all node ids to form a unique key.
    % We store the key as a 4-element row (pad with 0 for triangles).
    HF_key   = zeros(nHF, 4, 'int32');
    HF_verts = cell(nHF, 1);   % actual (unsorted) node list per half-face
    HF_cell  = zeros(nHF, 1, 'int32');
    HF_local = zeros(nHF, 1, 'int16');

    idx = 0;
    for ci = 1:nCells
        tp   = typs(ci);
        cnodes = cell_nodes{ci};
        lfd  = localFaceDef{tp};
        for lf = 1:numel(lfd)
            idx = idx + 1;
            local_ids = lfd{lf};
            vv = cnodes(local_ids);          % global node ids (1-based)
            sv = sort(double(vv));           % sorted key
            HF_key(idx, 1:numel(sv)) = int32(sv);
            HF_verts{idx} = vv;
            HF_cell(idx)  = ci;
            HF_local(idx) = lf;
        end
    end

    % Sort half-faces lexicographically to find matching pairs
    [HF_key_s, ia] = sortrows(HF_key);
    HF_verts_s = HF_verts(ia);
    HF_cell_s  = HF_cell(ia);

    % Group into unique faces: consecutive equal key rows are the same face
    nFaces      = 0;
    face_verts_c = cell(nHF, 1);   % will store node list for each unique face
    face_cells   = zeros(nHF, 2, 'int32');
    i = 1;
    while i <= nHF
        j = i;
        while j < nHF && all(HF_key_s(j+1,:) == HF_key_s(i,:))
            j = j + 1;
        end
        nFaces = nFaces + 1;
        face_verts_c{nFaces} = HF_verts_s{i};   % keep original winding
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
    %% 6.  For each cell: which faces does it own? (cell -> face map)
    %% ------------------------------------------------------------------ %%
    cell_faces_list = cell(nCells, 1);
    for f = 1:nFaces
        c1 = face_cells(f,1);
        c2 = face_cells(f,2);
        if c1 > 0
            cell_faces_list{c1} = [cell_faces_list{c1}, f];
        end
        if c2 > 0
            cell_faces_list{c2} = [cell_faces_list{c2}, f];
        end
    end

    %% ------------------------------------------------------------------ %%
    %% 7.  Compute face geometry: centroid, normal (unit), area
    %%
    %%  For polygonal faces (triangles or quads) we use the cross-product
    %%  decomposition: split into triangles from the first vertex.
    %% ------------------------------------------------------------------ %%
    face_centroid = zeros(nFaces, 3);
    face_normal   = zeros(nFaces, 3);
    face_area     = zeros(nFaces, 1);

    for f = 1:nFaces
        vids = double(face_verts_c{f});
        nv   = numel(vids);
        pts  = V3(vids, :);                  % nv x 3

        % Centroid = average of vertices
        face_centroid(f,:) = mean(pts, 1);

        % Area-weighted normal: sum of triangle contributions
        cr_sum = [0 0 0];
        p0 = pts(1,:);
        for kv = 2:nv-1
            cr_sum = cr_sum + cross(pts(kv,:) - p0, pts(kv+1,:) - p0);
        end
        area2 = norm(cr_sum);
        face_area(f) = area2 / 2;
        if area2 > 0
            face_normal(f,:) = cr_sum / area2;
        end
    end

    %% ------------------------------------------------------------------ %%
    %% 8.  Compute cell geometry: centroid and volume
    %%
    %%  General approach: decompose each cell into tetrahedra from its
    %%  centroid and accumulate signed volumes.  Works for any convex cell.
    %% ------------------------------------------------------------------ %%
    cell_centroid = zeros(nCells, 3);
    cell_volume   = zeros(nCells, 1);

    for ci = 1:nCells
        cnodes = double(cell_nodes{ci});
        % Approximate centroid as mean of vertices, then refine via tet decomp
        cpts   = V3(cnodes, :);
        xc_approx = mean(cpts, 1);

        tp  = typs(ci);
        lfd = localFaceDef{tp};
        vol  = 0;
        xc_w = [0 0 0];
        for lf = 1:numel(lfd)
            fvids = cnodes(lfd{lf});
            nfv   = numel(fvids);
            fpts  = V3(fvids, :);
            % Decompose face into triangles from first face vertex
            for kv = 2:nfv-1
                a  = fpts(1,:)   - xc_approx;
                b  = fpts(kv,:)  - xc_approx;
                c  = fpts(kv+1,:)- xc_approx;
                sv = dot(a, cross(b, c)) / 6;
                vol  = vol  + sv;
                % Tet centroid contribution (weight by signed volume)
                xc_w = xc_w + sv * (xc_approx + fpts(1,:) + fpts(kv,:) + fpts(kv+1,:)) / 4;
            end
        end
        cell_volume(ci)   = abs(vol);
        if abs(vol) > 0
            cell_centroid(ci,:) = xc_w / vol;
        else
            cell_centroid(ci,:) = xc_approx;
        end
    end

    %% ------------------------------------------------------------------ %%
    %% 9.  Populate face_struct  (matches MRSTGridConvert output exactly)
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
        face_struct(f).verts  = int32(face_verts_c{f}(:));  % column vector
        face_struct(f).center = face_centroid(f,:)';
        face_struct(f).normal = face_normal(f,:)';
        face_struct(f).area   = face_area(f);
    end

    %% ------------------------------------------------------------------ %%
    %% 10. Populate cell_struct  (matches MRSTGridConvert output exactly)
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

        % cells3D: unique vertex ids for this cell
        vl = double(cell_nodes{ci});
        cells3D{ci} = unique(vl(:), 'stable')';
    end

    fprintf('vtkLoader: done.\n');
    fprintf('  %d vertices\n', size(V3,1));
    fprintf('  %d cells\n', nCells);
    fprintf('  %d faces\n', nFaces);
end

%% ======================================================================
%% LOCAL HELPERS
%% ======================================================================

function val = xmlAttr(tag, attr)
% Return the string value of an XML attribute, or '' if not found.
    tok = regexpi(tag, [attr '\s*=\s*[''"]([^''"]*)[''"]'], 'tokens', 'once');
    if isempty(tok)
        val = '';
    else
        val = tok{1};
    end
end

function [da_tags, da_starts, da_ends] = findDataArrays(xml_str)
% Find all <DataArray ...>...</DataArray> blocks
    starts = regexp(xml_str, '<DataArray\s', 'start');
    da_tags   = cell(numel(starts),1);
    da_starts = zeros(numel(starts),1);
    da_ends   = zeros(numel(starts),1);
    for k = 1:numel(starts)
        % find end of opening tag
        tag_end = find(xml_str(starts(k):end) == '>', 1) + starts(k) - 1;
        da_tags{k}   = xml_str(starts(k):tag_end);
        da_starts(k) = tag_end + 1;
        close_pos    = strfind(xml_str(tag_end:end), '</DataArray>');
        da_ends(k)   = tag_end + close_pos(1) - 2;
    end
end

function raw_bytes = getDataByName(da_tags, da_starts, da_ends, xml_str, name)
% Extract and decode the DataArray with the given Name attribute.
% Returns raw decoded bytes (uint8), or [] if not found.
    raw_bytes = [];
    for k = 1:numel(da_tags)
        tok = regexpi(da_tags{k}, 'Name\s*=\s*[''"]([^''"]*)[''""]', 'tokens', 'once');
        if isempty(tok), continue; end
        if ~strcmpi(tok{1}, name), continue; end
        b64 = strtrim(xml_str(da_starts(k):da_ends(k)));
        % Strip any trailing XML child elements (e.g. <InformationKey>...)
        % Base64 data only contains A-Z, a-z, 0-9, +, /, =, and whitespace.
        % Anything starting with '<' is XML markup — truncate there.
        lt_pos = find(b64 == '<', 1);
        if ~isempty(lt_pos)
            b64 = strtrim(b64(1:lt_pos-1));
        end
        if isempty(b64), continue; end
        raw_bytes = decodeVTKBinary(b64);
        return;
    end
end

function data = decodeVTKBinary(b64str)
% Decode VTK binary block: base64 -> zlib decompress
%
% VTK with vtkZLibDataCompressor (header_type=UInt32) layout:
%   Base64( header ) ++ Base64( block_0 ) ++ ... ++ Base64( block_N-1 )
%
% Header (all uint32, little-endian):
%   [nblocks, uncompressed_block_size, last_partial_size,
%    compressed_size_0, ..., compressed_size_{N-1}]
%
% The header and each data block are encoded as *separate* base64 streams
% concatenated without any separator.  We identify boundaries by:
%   1. Decoding the first 4 bytes (= first 8 base64 chars, rounded to
%      multiple-of-4) to get nblocks.
%   2. The full header is (3 + nblocks)*4 bytes => known number of b64 chars.
%   3. Decode full header, read comp_sizes.
%   4. Decode rest as one flat byte array, slice by comp_sizes, decompress.

    % Strip all whitespace
    b64str = b64str(b64str ~= sprintf('\n') & b64str ~= sprintf('\r') & ...
                    b64str ~= ' '           & b64str ~= sprintf('\t'));

    if isempty(b64str)
        error('vtkLoader:decodeVTKBinary', 'Empty base64 string');
    end

    % --- Step 1: read nblocks from the first 4 bytes of the header.
    % 4 bytes needs ceil(4*4/3) = 6 b64 chars, padded to 8.
    hdr1 = b64decode_raw(b64str(1 : min(8, end)));
    nblocks = double(typecast(hdr1(1:4), 'uint32'));

    % --- Step 2: full header size and its base64 length.
    hdr_bytes_needed = (3 + nblocks) * 4;           % bytes
    hdr_b64_len = ceil(hdr_bytes_needed / 3) * 4;   % base64 chars (no padding yet)

    hdr_b64 = b64str(1 : hdr_b64_len);
    hdr_all = b64decode_raw(hdr_b64);
    hdr_all = hdr_all(:)';  % row vector

    comp_sizes = zeros(1, nblocks, 'double');
    for b = 1:nblocks
        off = 12 + (b-1)*4 + 1;
        comp_sizes(b) = double(typecast(hdr_all(off : off+3), 'uint32'));
    end

    % --- Step 3: decode the rest as flat bytes.
    rest_b64 = b64str(hdr_b64_len + 1 : end);
    all_data = b64decode_raw(rest_b64);
    all_data = all_data(:)';  % row vector

    % --- Step 4: slice and decompress each block.
    pos  = 1;
    data = uint8([]);
    for b = 1:nblocks
        cs = comp_sizes(b);
        if pos + cs - 1 > numel(all_data)
            error('vtkLoader:decodeVTKBinary', ...
                'Block %d: need %d bytes at pos %d but data has only %d bytes', ...
                b, cs, pos, numel(all_data));
        end
        chunk        = all_data(pos : pos+cs-1);
        decompressed = zlib_decompress(chunk(:));
        data         = [data, decompressed(:)']; %#ok<AGROW>
        pos          = pos + cs;
    end
    data = data(:);  % column vector for typecast callers
end

function out = b64decode_raw(str)
% Decode a base64 string (padding optional) to uint8 row vector.
% Does NOT require the string to be padded to a multiple of 4.
    table = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
    % Strip padding
    str = str(str ~= '=');
    n   = numel(str);
    if n == 0
        out = uint8([]);
        return;
    end
    % Map characters to 6-bit values
    vals = zeros(1, n, 'uint8');
    for k = 1:n
        idx = find(table == str(k), 1);
        if isempty(idx)
            error('vtkLoader:b64decode_raw', 'Invalid base64 char: %s', str(k));
        end
        vals(k) = idx - 1;
    end
    % Decode groups
    nGroups = floor(n / 4);
    rem_    = mod(n, 4);
    nout    = nGroups * 3;
    if rem_ == 2,      nout = nout + 1;
    elseif rem_ == 3,  nout = nout + 2; end
    out = zeros(1, nout, 'uint8');
    oi  = 1;
    for g = 1:nGroups
        v = vals((g-1)*4+1 : g*4);
        bits = bitor(bitor(bitor( ...
            bitshift(uint32(v(1)), 18), bitshift(uint32(v(2)), 12)), ...
            bitshift(uint32(v(3)),  6)), uint32(v(4)));
        out(oi)   = uint8(bitshift(bits, -16));
        out(oi+1) = uint8(bitand(bitshift(bits, -8), 255));
        out(oi+2) = uint8(bitand(bits, 255));
        oi = oi + 3;
    end
    % Remainder
    if rem_ >= 2
        v1 = vals(nGroups*4+1);
        v2 = vals(nGroups*4+2);
        out(oi) = uint8(bitor(bitshift(uint32(v1), 2), bitshift(uint32(v2), -4)));
        oi = oi + 1;
    end
    if rem_ == 3
        v2 = vals(nGroups*4+2);
        v3 = vals(nGroups*4+3);
        out(oi) = uint8(bitand(bitor(bitshift(uint32(v2), 4), bitshift(uint32(v3), -2)), 255));
    end
end

function tokens = splitBase64Tokens(str)
% Split a concatenated base64 string at padding boundaries.
% VTK writes each block as a separate base64 stream that ends with '='.
    % Strip all whitespace first
    str = str(str ~= sprintf('\n') & str ~= sprintf('\r') & ...
               str ~= ' '          & str ~= sprintf('\t'));

    tokens = {};
    start  = 1;
    n      = numel(str);
    i      = 1;
    while i <= n
        if str(i) == '='
            % Consume run of '=' (at most 2 for valid base64)
            j = i;
            while j <= n && str(j) == '='
                j = j + 1;
            end
            tok = str(start:j-1);
            r = mod(numel(tok), 4);
            if r ~= 0
                tok = [tok, repmat('=', 1, 4-r)]; %#ok<AGROW>
            end
            tokens{end+1} = tok; %#ok<AGROW>
            start = j;
            i     = j;
        else
            i = i + 1;
        end
    end
    % Any remaining characters after the last '='
    if start <= n
        tok = str(start:n);
        r   = mod(numel(tok), 4);
        if r ~= 0
            tok = [tok, repmat('=', 1, 4-r)];
        end
        tokens{end+1} = tok; %#ok<AGROW>
    end
end

function out = base64decode_chunk(str)
% Decode a single, complete base64 stream (may end with 0,1,2 '=').
    table = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
    % Pad to multiple of 4 if needed
    r = mod(numel(str), 4);
    if r ~= 0
        str = [str, repmat('=', 1, 4-r)];
    end
    n    = numel(str);
    vals = zeros(1, n, 'uint8');
    for k = 1:n
        c = str(k);
        if c == '='
            vals(k) = 0;
        else
            idx = find(table == c, 1);
            if isempty(idx)
                error('vtkLoader:base64decode', 'Unexpected character ''%s'' in base64 stream', c);
            end
            vals(k) = idx - 1;
        end
    end
    nGroups = n / 4;
    out = zeros(1, nGroups * 3, 'uint8');
    for g = 1:nGroups
        v = vals((g-1)*4+1 : g*4);
        bits = bitor(bitor(bitor( ...
            bitshift(uint32(v(1)), 18), ...
            bitshift(uint32(v(2)), 12)), ...
            bitshift(uint32(v(3)),  6)), ...
            uint32(v(4)));
        out((g-1)*3+1) = uint8(bitshift(bits, -16));
        out((g-1)*3+2) = uint8(bitand(bitshift(bits, -8), 255));
        out((g-1)*3+3) = uint8(bitand(bits, 255));
    end
    % Remove padding bytes
    if numel(str) >= 2
        npad = sum(str(end-1:end) == '=');
    elseif numel(str) == 1
        npad = (str(end) == '=');
    else
        npad = 0;
    end
    out  = out(1 : end - npad);
end

function out = zlib_decompress(data)
% Decompress a zlib-compressed byte vector using Java.
% Tries standard zlib first; falls back to raw deflate (nowrap=true).
% Always returns a row vector.
    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier();
    try
        a = java.io.ByteArrayInputStream(data);
        b = java.util.zip.InflaterInputStream(a);
        c = java.io.ByteArrayOutputStream();
        isc.copyStream(b, c);
        out = typecast(c.toByteArray(), 'uint8')';
        out = out(:)';  % ensure row
    catch
        % Try raw deflate (no zlib wrapper)
        try
            inf = java.util.zip.Inflater(true);   % nowrap = true
            inf.setInput(data);
            buf = java.io.ByteArrayOutputStream();
            tmp = zeros(1, 65536, 'int8');
            while ~inf.finished()
                n = inf.inflate(tmp);
                if n > 0
                    buf.write(tmp, 0, n);
                end
                if n == 0 && ~inf.finished()
                    break;
                end
            end
            inf.end();
            out = typecast(buf.toByteArray(), 'uint8')';
            out = out(:)';  % ensure row
        catch ME2
            error('vtkLoader:zlib_decompress', ...
                'Could not decompress block (tried zlib and raw deflate): %s', ME2.message);
        end
    end
end

