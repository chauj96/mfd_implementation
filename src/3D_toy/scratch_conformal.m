function G = scratch_conformal()
mrstModule add upr
%% set boundary
tx = 1;
ty = 1;
tz = 0.4;
bdr   = [ 0, 0, 0;  ...
          tx, 0, 0;  ...
          tx, ty, 0;  ...
          0, ty, 0;  ...
          0, 0, tz;  ...
          tx, 0, tz;  ...
          tx, ty, tz;  ...
          0, ty, tz];

%% Set gridding parameters
fGs = tx/16;
rho = @(p) fGs*(1+0*p(:,1));

%% Triangulate the two curved surfaces
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
rectangle1 = [min(bdr(:, [1,3])); max(bdr(:,[1,3]))];
fixedPts = [rectangle1; rectangle1(swap1)'; rectangle1(swap2)'];
hd = @(p) rho(p)/fGs;
fd = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3), rectangle1(4));

[Pts,t] = distmesh2d(fd, hd, fGs, rectangle1, fixedPts, false);

t1.ConnectivityList = t;
t2.ConnectivityList = t;

% Surface geometry parameters — reduce these to make surfaces flatter:
%   zTilt    : vertical tilt across y  (0 = horizontal plane)
%   yBend    : parabolic bend in y     (0 = flat in y)
%   sinAmp   : sinusoidal waviness     (0 = no waviness)
zTilt  = 0.15;   % was 0.15
yBend  = 0.0;   % was 0.08
sinAmp = 0.0;   % was 0.04

faultHeight1z = @(p)   ty/6*ones(size(p,1),1) + zTilt*p(:,2);
faultHeight2z = @(p) 5*ty/6*ones(size(p,1),1) - zTilt*p(:,2);

t1.Points = [Pts(:,1), faultHeight1z(Pts), Pts(:,2)];
t2.Points = [Pts(:,1), faultHeight2z(Pts), Pts(:,2)];

faultHeight1x = @(p) p(:,2) + yBend*(p(:,1)-0.75).^2 + sinAmp*sin(2*pi/tx*p(:,1));
faultHeight2x = @(p) p(:,2) + yBend*(p(:,1)-0.75).^2;

t1.Points(:,2) = faultHeight1x(t1.Points);
t2.Points(:,2) = faultHeight2x(t2.Points);

%% Generate conformal surface sites using the UPR pipeline
% surfaceSites3D places equidistant site pairs on both sides of each
% triangulated surface so that the Voronoi faces align with the surfaces.
%
% Radius strategy: per-vertex R = max_incident_circumradius / gamma.
%
% ballInt places each site pair at the circumcenter of its triangle,
% offset by Z = sqrt(R^2 - Rcirc^2) along the surface normal.
% We need R > Rcirc for every incident triangle, so we set
%   R = max(Rcirc of incident triangles) / gamma,   gamma < 1
% which guarantees Z is real and positive.  Smaller gamma -> larger R
% -> larger Z (wider site separation) -> better-shaped 3D Voronoi cells.
% gamma = R / max_incident_Rcirc.  Must be > 1 so R > Rcirc and Z is real.
% gamma=0.9 (i.e. R = Rcirc/0.9) gave the best inR/outR in empirical tests.
% Do NOT use the mean Rcirc to derive gamma: Gs varies per-triangle, so a
% mean-based formula produces wildly inconsistent Gs values across the mesh.
gamma = 0.9;

R1 = circumradiusPerVertex(t1.Points, t1.ConnectivityList, gamma);
R2 = circumradiusPerVertex(t2.Points, t2.ConnectivityList, gamma);

fprintf('\n=== Circumradius-based ball radius (Surface 1) ===\n');
fprintf('  R min/mean/max: %.5f / %.5f / %.5f\n', min(R1), mean(R1), max(R1));
fprintf('=== Circumradius-based ball radius (Surface 2) ===\n');
fprintf('  R min/mean/max: %.5f / %.5f / %.5f\n', min(R2), mean(R2), max(R2));

rho1 = @(p) R1;
rho2 = @(p) R2;

F = surfaceSites3D({t1, t2}, {rho1, rho2});

%% Hard filter: remove any complex, NaN, or Inf sites produced by ballInt
isValid = all(isreal(F.f.pts), 2) & all(isfinite(F.f.pts), 2);
nBad = sum(~isValid);
if nBad > 0
    fprintf('\n  WARNING: removed %d complex/NaN surface sites from ballInt.\n', nBad);
end
F.f.pts = F.f.pts(isValid, :);
F.f.Gs  = F.f.Gs(isValid);
F.f.pri = F.f.pri(isValid);
F.f.l   = F.f.l(isValid);

fprintf('  Surface site pairs: %d\n', size(F.f.pts,1));
fprintf('  Site-pair separation Gs  min/mean/max: %.4f / %.4f / %.4f\n', ...
        min(F.f.Gs), mean(F.f.Gs), max(F.f.Gs));
fprintf('  Gs / fGs ratio           min/mean/max: %.3f / %.3f / %.3f\n', ...
        min(F.f.Gs)/fGs, mean(F.f.Gs)/fGs, max(F.f.Gs)/fGs);

%% Create background reservoir sites on a regular grid
dt = fGs; % use the same mesh size
xmax = max(bdr(:,1)) - dt/2;  xmin = min(bdr(:,1)) + dt/2;
ymax = max(bdr(:,2)) - dt/2;  ymin = min(bdr(:,2)) + dt/2;
zmax = max(bdr(:,3)) - dt/2;  zmin = min(bdr(:,3)) + dt/2;

xr = xmin:dt:xmax;  yr = ymin:dt:ymax;  zr = zmin:dt:zmax;
[X,Y,Z] = ndgrid(xr, yr, zr);
rSites = [X(:), Y(:), Z(:)];

% Jitter background sites slightly to break exact equidistance configurations.
rng(42);
rSites = rSites + 0.01*dt*(2*rand(size(rSites))-1);

% Step 1 – necessary condition: remove sites inside circumscribed vertex balls.
rSites = surfaceSufCond3D(rSites, F.c.CC, F.c.R);

% Step 2 – sufficient condition: remove reservoir sites within F.f.Gs of
% each surface site.  F.f.Gs is the actual site-pair separation for that
% triangle; using it as the exclusion radius ensures no reservoir site can
% be closer to a surface site than the surface site is to its paired copy,
% which is the exact condition that prevents flat bisecting Voronoi faces.
rSites = surfaceSufCond3D(rSites, F.f.pts, F.f.Gs);

fprintf('  Reservoir sites after exclusion: %d\n', size(rSites,1));

%% Merge surface sites and reservoir sites, then build grid
sites = [F.f.pts; rSites];

% Clip all sites strictly inside bdr so mirroredPebi3D does not silently
% drop surface sites that landed exactly on or outside the boundary.
eps_clip = 1e-8;
sites(:,1) = min(max(sites(:,1), min(bdr(:,1)) + eps_clip), max(bdr(:,1)) - eps_clip);
sites(:,2) = min(max(sites(:,2), min(bdr(:,2)) + eps_clip), max(bdr(:,2)) - eps_clip);
sites(:,3) = min(max(sites(:,3), min(bdr(:,3)) + eps_clip), max(bdr(:,3)) - eps_clip);

fprintf('Total number of sites = %d\n', size(sites,1));

G = mirroredPebi3D(sites, bdr);
G = computeGeometry(G);

% Compute robust geometric centroids for all cells (replaces G.cells.centroids
% which holds Voronoi generators for surface cells after mirroredPebi3D).
geom_cc  = zeros(G.cells.num, 3);
geom_vol = zeros(G.cells.num, 1);
for ci = 1:G.cells.num
    [geom_cc(ci,:), geom_vol(ci)] = computeRobustCellCentroid(G, ci);
end

surf_sites = F.f.pts;
nSurf      = size(surf_sites, 1);
nCells_g   = G.cells.num;


G.cells.centroids = geom_cc;


fprintf('Mesh built:\n');
fprintf('  #cells = %d\n', G.cells.num);
fprintf('  #faces = %d\n', G.faces.num);
fprintf('  #nodes = %d\n', G.nodes.num);

%% Classify cells by position relative to the two surfaces
c  = G.cells.centroids;
cl = c(:,2) > faultHeight2z(c(:,[1,3])) + faultHeight2x(c(:,[1,3])) - c(:,3);
cr = c(:,2) < faultHeight1z(c(:,[1,3])) + faultHeight1x(c(:,[1,3])) - c(:,3);
cm = ~cl & ~cr;

%% Plot
figure(); hold on
color = lines(7);

plotGrid(G, cl, 'FaceColor', color(1,:), 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.15);
plotGrid(G, cm, 'FaceColor', color(2,:), 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.15);
plotGrid(G, cr, 'FaceColor', color(3,:), 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.15);
plotGrid(G,      'FaceColor', 'none',    'EdgeColor', [0 0 0], 'EdgeAlpha', 0.2);

trisurf(t1.ConnectivityList, ...
    t1.Points(:,1), t1.Points(:,2), t1.Points(:,3), ...
    'FaceColor', color(1,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
trisurf(t2.ConnectivityList, ...
    t2.Points(:,1), t2.Points(:,2), t2.Points(:,3), ...
    'FaceColor', color(2,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');

view(120,30);
light('Position',[5 10 -1],'Style','infinite')
light('Position',[-1 -1 2],'Style','local')
axis equal off
set(gca,'zdir','normal')
title('Conformal PEBI mesh with two curved surfaces')

%% Euler characteristic  χ = V - E + F - C
% For a 3-D cell complex: χ = #vertices - #edges + #faces - #cells
% MRST stores vertices (nodes) and faces explicitly; edges must be derived
% from the face-node connectivity list.
V = G.nodes.num;
F_count = G.faces.num;
C = G.cells.num;

% Build edge list: every face is a polygon; collect all consecutive
% node pairs (including the closing edge) and keep unique undirected ones.
faceNodes  = G.faces.nodes;          % concatenated node indices
faceNCount = diff(G.faces.nodePos);  % number of nodes per face

% Build start/end node pairs for every edge of every face
startNodes = faceNodes;
% For each face, shift the node list by one to get the next node;
% the last node of each face connects back to its first.
endNodes = zeros(size(faceNodes), 'like', faceNodes);
for fIdx = 1:F_count
    nStart = G.faces.nodePos(fIdx);
    nEnd   = G.faces.nodePos(fIdx+1) - 1;
    idx    = nStart:nEnd;
    endNodes(idx) = faceNodes([idx(2:end), idx(1)]);
end

% Make edges undirected (smaller index first) and deduplicate
edges = sort([startNodes, endNodes], 2);
edges = unique(edges, 'rows');
E = size(edges, 1);

chi = V - E + F_count - C;

fprintf('\nEuler characteristic:\n');
fprintf('  V (vertices) = %d\n', V);
fprintf('  E (edges)    = %d\n', E);
fprintf('  F (faces)    = %d\n', F_count);
fprintf('  C (cells)    = %d\n', C);
fprintf('  chi = V - E + F - C = %d - %d + %d - %d = %d\n', ...
        V, E, F_count, C, chi);

%% --- Volume consistency check ---
checkVolumeConsistency(G);

end


%% ---------------------------------------------------------------
function checkVolumeConsistency(G)
nCells = G.cells.num;
fc     = G.faces.centroids;
fn     = G.faces.normals;
cc     = G.cells.centroids;
cv     = G.cells.volumes;
N      = G.faces.neighbors;

cFacePos = G.cells.facePos;
cFaces   = G.cells.faces(:,1);

inexact     = 0;
max_rel_err = 0;
fail_list   = [];

for ci = 1:nCells
    idx  = cFacePos(ci) : cFacePos(ci+1)-1;
    fids = cFaces(idx);
    xc   = cc(ci,:)';
    Vc   = cv(ci);

    T = zeros(3,3);
    for k = 1:numel(fids)
        f   = fids(k);
        sgn = 2*(N(f,1) == ci) - 1;   % +1 or -1
        n_out = sgn * fn(f,:)';
        xf    = fc(f,:)';
        T     = T + n_out * (xf - xc)';
    end

    % Also recompute volume from divergence theorem (trace/3)
    V_div = trace(T) / 3;

    rel_err     = norm(T - Vc * eye(3), 'fro') / Vc;
    rel_vol_err = abs(V_div - Vc) / Vc;
    max_rel_err = max(max_rel_err, rel_err);
    if rel_err > 1e-10
        inexact = inexact + 1;
        fail_list(end+1,:) = [ci, rel_err, Vc, V_div, rel_vol_err]; %#ok<AGROW>
    end
end

fprintf('\n=== scratch_conformal:: MFD Volume Consistency Check ===\n');
fprintf('  Cells checked           : %d\n', nCells);
fprintf('  Max relative error      : %.3e\n', max_rel_err);
if inexact == 0
    fprintf('  All cells PASS (tol 1e-10).\n');
else
    fprintf('  FAILED cells (rel_err > 1e-10): %d / %d\n', inexact, nCells);

    % Diagnose: are failures concentrated near face centroids?
    % Recompute volume via div theorem and compare with stored volume
    vol_err = fail_list(:,5);
    fprintf('\n  Volume discrepancy (|V_div - V_stored|/V_stored) in failed cells:\n');
    fprintf('    min/mean/max: %.3e / %.3e / %.3e\n', ...
            min(vol_err), mean(vol_err), max(vol_err));
    fprintf('  => If vol discrepancy is large, centroid correction broke\n');
    fprintf('     the centroid/volume consistency assumed by computeGeometry.\n');

    % Show worst 5
    [~, ord] = sort(fail_list(:,2), 'descend');
    fprintf('\n  Worst 5 cells (ci, rel_T_err, V_stored, V_div):\n');
    for k = 1:min(5, size(fail_list,1))
        r = fail_list(ord(k),:);
        fprintf('    cell %4d  rel_err=%.3e  V_stored=%.4e  V_div=%.4e\n', ...
                r(1), r(2), r(3), r(4));
    end
end
end


function [centroid, volume] = computeRobustCellCentroid(G, ci)
% computeRobustCellCentroid  Compute cell centroid and volume via the
%   divergence-theorem / signed-tetrahedra method.
%
%   For each polygonal face, fan-triangulate from the face centroid,
%   form a tetrahedron with an arbitrary reference point, accumulate
%   signed volumes and centroid contributions.
%
%   This is numerically stable even for:
%     - non-convex cells
%     - boundary cells far from the origin
%     - thin/sliver conformal surface cells
%
%   Inputs
%     G   : MRST grid (computeGeometry already called)
%     ci  : cell index
%
%   Outputs
%     centroid : 1x3 centroid (geometric, not Voronoi generator)
%     volume   : signed volume (positive for outward-consistent normals)

    % ---- gather face indices for this cell --------------------------
    faceIdx = G.cells.faces(G.cells.facePos(ci):G.cells.facePos(ci+1)-1, 1);
    nFaces  = numel(faceIdx);

    % ---- reference point: mean of face centroids (avoids origin bias) --
    ref = mean(G.faces.centroids(faceIdx, :), 1);   % 1x3

    v_total    = 0.0;
    c_weighted = [0.0, 0.0, 0.0];

    for k = 1:nFaces
        f = faceIdx(k);

        % Outward sign: MRST normals point from neighbor(1)->neighbor(2)
        sgn = 2 * (G.faces.neighbors(f, 1) == ci) - 1;   % +1 or -1

        % Nodes of this face
        nStart = G.faces.nodePos(f);
        nEnd   = G.faces.nodePos(f+1) - 1;
        fnodes = G.faces.nodes(nStart:nEnd);   % polygon node indices
        nv     = numel(fnodes);

        % Face centroid as fan apex (handles non-planar/non-convex faces)
        fc = mean(G.nodes.coords(fnodes, :), 1);   % 1x3

        % Fan-triangulate the polygon from fc
        for j = 1:nv
            % Triangle vertices (shifted to ref for stability)
            a = G.nodes.coords(fnodes(j),              :) - ref;   % 1x3
            b = G.nodes.coords(fnodes(mod(j,nv)+1),    :) - ref;   % 1x3
            c = fc - ref;                                           % 1x3

            % Ensure outward orientation
            if sgn < 0
                [a, b] = deal(b, a);   % swap to flip triangle orientation
            end

            % Signed tetrahedron volume: (1/6) * dot(a, cross(b,c))
            v_i = (1.0/6.0) * dot(a, cross(b, c));

            % Centroid of this tetrahedron: (ref + a + b + c) / 4
            % But since a,b,c are already shifted, tet centroid = (a+b+c)/4
            v_total    = v_total    + v_i;
            c_weighted = c_weighted + v_i * (a + b + c) / 4.0;
        end
    end

    volume   = v_total;
    centroid = ref + c_weighted / v_total;   % shift back to world coords
end


function [face_centroid, face_area, face_normal] = computeRobustFaceCentroid(G, f)
% computeRobustFaceCentroid  Robust face centroid, area, and unit normal
%   via fan-triangulation from the vertex mean.
%
%   Method: fan-triangulate polygon from vertex mean (ref), accumulate
%   area-weighted sub-triangle centroids. Consistent with the signed-
%   tetrahedra cell centroid method.
%
%   Inputs
%     G : MRST grid (computeGeometry already called)
%     f : face index (1-based)
%
%   Outputs
%     face_centroid : 1x3  area-weighted centroid
%     face_area     : scalar >= 0
%     face_normal   : 1x3  unit normal, sign consistent with MRST
%                     (points from neighbor(1) toward neighbor(2))

    % --- gather polygon nodes -------------------------------------------
    nStart = G.faces.nodePos(f);
    nEnd   = G.faces.nodePos(f+1) - 1;
    fnodes = G.faces.nodes(nStart:nEnd);
    nv     = numel(fnodes);
    coords = G.nodes.coords(fnodes, :);   % nv x 3

    % --- reference point: vertex mean (removes origin bias) -------------
    ref = mean(coords, 1);                % 1x3

    % --- fan-triangulation from ref -------------------------------------
    area_vec_sum = [0.0, 0.0, 0.0];      % signed area-weighted normal
    c_weighted   = [0.0, 0.0, 0.0];      % area-weighted centroid accumulator
    area_total   = 0.0;                   % sum of unsigned triangle areas

    for j = 1:nv
        a = coords(j,           :) - ref;              % 1x3
        b = coords(mod(j,nv)+1, :) - ref;              % 1x3
        % apex is ref itself => c = ref - ref = [0,0,0]

        cr       = cross(a, b);                         % signed area vector * 2
        tri_area = 0.5 * norm(cr);

        if tri_area < eps
            continue                                    % skip degenerate sliver
        end

        area_vec_sum = area_vec_sum + 0.5 * cr;         % accumulate signed
        area_total   = area_total   + tri_area;

        % Sub-triangle centroid in shifted coords: (a + b + 0) / 3
        c_weighted = c_weighted + tri_area * (a + b) / 3.0;
    end

    % --- face area & centroid -------------------------------------------
    face_area = norm(area_vec_sum);                     % signed magnitude

    if area_total > eps
        face_centroid = ref + c_weighted / area_total;
    else
        % Degenerate face: fall back to vertex mean
        face_centroid = ref;
        warning('computeRobustFaceCentroid: degenerate face %d, using vertex mean', f);
    end

    % --- unit normal, sign consistent with MRST -------------------------
    % MRST convention: normals(f) points from neighbor(1) -> neighbor(2).
    % area_vec_sum direction depends on node winding; align with MRST.
    if face_area > eps
        n_candidate = area_vec_sum / face_area;

        % Validate against MRST's stored normal (already area-weighted)
        n_mrst = G.faces.normals(f, :);
        n_mrst_unit = n_mrst / max(norm(n_mrst), eps);

        if dot(n_candidate, n_mrst_unit) < 0
            n_candidate = -n_candidate;               % flip to match MRST
        end
        face_normal = n_candidate;
    else
        % Degenerate: use MRST stored normal as fallback
        n_mrst = G.faces.normals(f, :);
        n_norm = norm(n_mrst);
        face_normal = n_mrst / max(n_norm, eps);
        warning('computeRobustFaceCentroid: near-zero area face %d', f);
    end
end


%% ---------------------------------------------------------------
function R = circumradiusPerVertex(pts, conn, gamma)
% Per-vertex ball radius = max_incident_circumradius / gamma.
% ballInt computes Z = sqrt(R^2 - Rcirc^2); we need R > Rcirc for every
% incident triangle, so R = max(Rcirc) / gamma with gamma < 1.
nPts = size(pts, 1);
p1 = pts(conn(:,1),:);  p2 = pts(conn(:,2),:);  p3 = pts(conn(:,3),:);
L1 = sqrt(sum((p2-p1).^2,2));
L2 = sqrt(sum((p3-p2).^2,2));
L3 = sqrt(sum((p1-p3).^2,2));
s     = (L1+L2+L3)/2;
area  = sqrt(max(s.*(s-L1).*(s-L2).*(s-L3), 0));
Rcirc = (L1.*L2.*L3) ./ max(4*area, eps);
verts = [conn(:,1); conn(:,2); conn(:,3)];
rvals = [Rcirc;     Rcirc;     Rcirc    ];
R = accumarray(verts, rvals, [nPts,1], @max);
R = R / gamma;
end




