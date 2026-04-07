function G = scratch()
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
fGs = tx/20;
rho = @(p) fGs*(1+0*p(:,1));

%% Triangulate the two curved surfaces
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
rectangle1 = [min(bdr(:, [1,3])); max(bdr(:,[1,3]))];
fixedPts = [rectangle1; rectangle1(swap1)'; rectangle1(swap2)'];
hd = @(p) rho(p)/fGs;
fd = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3), rectangle1(4));

[Pts,t] = distmesh2d(fd, hd, fGs, rectangle1, fixedPts, false);

%% --- Diagnostic: 2D parameter-space triangle quality ---
fprintf('\n=== Surface triangulation quality (2D parameter space) ===\n');
triQualityReport(Pts, t, '2D param space');

t1.ConnectivityList = t;
t2.ConnectivityList = t;

faultHeight1z = @(p)   ty/6*ones(size(p,1),1) + 0.4*p(:,2);
faultHeight2z = @(p) 5*ty/6*ones(size(p,1),1) - 0.4*p(:,2);

t1.Points = [Pts(:,1), faultHeight1z(Pts), Pts(:,2)];
t2.Points = [Pts(:,1), faultHeight2z(Pts), Pts(:,2)];

faultHeight1x = @(p) p(:,2) + 0.2*(p(:,1)-0.75).^2 + 0.1*sin(2*pi/tx*p(:,1));
faultHeight2x = @(p) p(:,2) + 0.2*(p(:,1)-0.75).^2;

t1.Points(:,2) = faultHeight1x(t1.Points);
t2.Points(:,2) = faultHeight2x(t2.Points);

%% --- Diagnostic: 3D surface quality for both surfaces ---
fprintf('\n=== Surface 1 triangle quality (3D) ===\n');
triQualityReport(t1.Points, t1.ConnectivityList, 'Surface 1 (3D)');

fprintf('\n=== Surface 2 triangle quality (3D) ===\n');
triQualityReport(t2.Points, t2.ConnectivityList, 'Surface 2 (3D)');

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
gamma = 0.6;   % 0 < gamma < 1; smaller -> larger Z -> rounder surface cells

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

%% Expand boundary slightly for robustness
bdr = 1.01 * bdr;

%% Create background reservoir sites on a regular grid
dt = fGs; % use the same mesh size
xmax = max(bdr(:,1)) - dt/2;  xmin = min(bdr(:,1)) + dt/2;
ymax = max(bdr(:,2)) - dt/2;  ymin = min(bdr(:,2)) + dt/2;
zmax = max(bdr(:,3)) - dt/2;  zmin = min(bdr(:,3)) + dt/2;

xr = xmin:dt:xmax;  yr = ymin:dt:ymax;  zr = zmin:dt:zmax;
[X,Y,Z] = ndgrid(xr, yr, zr);
rSites = [X(:), Y(:), Z(:)];

% Jitter background sites slightly to break exact equidistance
% configurations that produce zero-area Voronoi faces.
rng(42);
rSites = rSites + 0.02*dt*(2*rand(size(rSites))-1);

% Step 1: remove reservoir sites inside the circumscribed balls centred at
% each surface triangulation vertex (necessary conformity condition).
rSites = surfaceSufCond3D(rSites, F.c.CC, F.c.R);

% Step 2: remove reservoir sites within Gs/2 of any surface site.
% F.f.pts are the conformal site pairs; F.f.Gs is the separation of each
% pair.  A reservoir site closer than Gs/2 to either site of a pair would
% have that surface site as its Voronoi neighbour, producing a flat cell.
% We use surfaceSufCond3D with the surface sites as centres and Gs/2 as radii.
halfGs = F.f.Gs / 2;
rSites = surfaceSufCond3D(rSites, F.f.pts, halfGs);

fprintf('  Reservoir sites after exclusion: %d\n', size(rSites,1));

%% Merge surface sites and reservoir sites, then build grid
sites = [F.f.pts; rSites];

fprintf('Total number of sites = %d\n', size(sites,1));

G = mirroredPebi3D(sites, bdr);

%% Clean up degenerate geometry using removeShortEdges (public API)
%
% Sliver faces have area << mean but topologically valid nodes.
% Their shortest edge has length ≈ 2*area / longest_edge ≈ 2*area / fGs.
% Setting edgeTol just above that length causes removeShortEdges to merge
% the two nodes of that short edge, which collapses the sliver face to
% < 3 nodes → removeCollapsedFaces removes it → removeCollapsedCells
% merges the two formerly-separated cells.
%
% Threshold: 1e-4 * mean(face areas) — conservative enough to only touch
% genuinely degenerate faces, not valid small-but-real faces.
G_tmp   = computeGeometry(G);
areaTol = 1e-4 * mean(G_tmp.faces.areas);
edgeTol = 2 * areaTol / fGs;
clear G_tmp;

fprintf('\nCleaning degenerate geometry (edgeTol = %.2e) ...\n', edgeTol);
nC0 = G.cells.num; nF0 = G.faces.num; nN0 = G.nodes.num;
[G, ~] = removeShortEdges(G, edgeTol);
fprintf('  Removed: %d cells,  %d faces,  %d nodes\n', ...
    nC0-G.cells.num, nF0-G.faces.num, nN0-G.nodes.num);


G = computeGeometry(G);

fprintf('Mesh built:\n');
fprintf('  #cells = %d\n', G.cells.num);
fprintf('  #faces = %d\n', G.faces.num);
fprintf('  #nodes = %d\n', G.nodes.num);

%% --- Diagnostic: 3D Voronoi cell quality ---
cellQualityReport(G);

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

%% ---------------------------------------------------------------
function cellQualityReport(G)
% Report quality of 3D Voronoi cells: insphere/circumsphere ratio,
% volume distribution, and face area distribution.
cc = G.cells.centroids;
nC = G.cells.num;
fc       = G.faces.centroids;
cFacePos = G.cells.facePos;
inR  = zeros(nC,1);
outR = zeros(nC,1);
for k = 1:nC
    fi      = G.cells.faces(cFacePos(k):cFacePos(k+1)-1);
    d       = sqrt(sum((fc(fi,:) - cc(k,:)).^2, 2));
    inR(k)  = min(d);
    outR(k) = max(d);
end
ratio = inR ./ max(outR, eps);
vol   = G.cells.volumes;
area  = G.faces.areas;
fprintf('\n=== 3D Voronoi cell quality ===\n');
fprintf('  Cell volume  min/mean/max : %.2e / %.2e / %.2e\n', min(vol),mean(vol),max(vol));
fprintf('  Volume CoV (std/mean)     : %.3f\n', std(vol)/mean(vol));
fprintf('  inR/outR  min/mean/max    : %.3f / %.3f / %.3f\n', min(ratio),mean(ratio),max(ratio));
fprintf('  Frac. inR/outR < 0.10     : %.1f%%  (flat/sliver cells)\n', 100*mean(ratio<0.10));
fprintf('  Frac. inR/outR < 0.20     : %.1f%%\n', 100*mean(ratio<0.20));
fprintf('  Frac. vol < 0.01*mean_vol : %.1f%%  (tiny cells)\n', 100*mean(vol<0.01*mean(vol)));
fprintf('  Face area  min/mean/max   : %.2e / %.2e / %.2e\n', min(area),mean(area),max(area));
fprintf('  Frac. face area < 1e-10   : %.2f%%  (degenerate faces)\n', 100*mean(area<1e-10));
fprintf('  Frac. face area < 1e-6    : %.2f%%\n', 100*mean(area<1e-6));
end

%% ---------------------------------------------------------------
function triQualityReport(pts, conn, label)
% Triangle quality: edge lengths, aspect ratio, minimum angle.
nT = size(conn,1);
p1 = pts(conn(:,1),:);  p2 = pts(conn(:,2),:);  p3 = pts(conn(:,3),:);
L1 = sqrt(sum((p2-p1).^2,2));
L2 = sqrt(sum((p3-p2).^2,2));
L3 = sqrt(sum((p1-p3).^2,2));
allL = [L1;L2;L3];
s    = (L1+L2+L3)/2;
area = sqrt(max(s.*(s-L1).*(s-L2).*(s-L3), 0));
R_circ = (L1.*L2.*L3) ./ max(4*area, eps);
r_in   = area ./ max(s, eps);
AR     = R_circ ./ max(2*r_in, eps);
cosA = (L2.^2+L3.^2-L1.^2) ./ max(2*L2.*L3, eps);
cosB = (L1.^2+L3.^2-L2.^2) ./ max(2*L1.*L3, eps);
cosC = (L1.^2+L2.^2-L3.^2) ./ max(2*L1.*L2, eps);
minAng = min([acosd(min(max(cosA,-1),1)), ...
              acosd(min(max(cosB,-1),1)), ...
              acosd(min(max(cosC,-1),1))], [], 2);
nDegen = sum(area < 1e-14);
fprintf('[%s]\n', label);
fprintf('  #triangles          : %d  (degenerate: %d)\n', nT, nDegen);
fprintf('  Edge length  min/mean/max: %.4f / %.4f / %.4f\n', min(allL),mean(allL),max(allL));
fprintf('  Aspect ratio min/mean/max: %.3f / %.3f / %.3f\n', min(AR),mean(AR),max(AR));
fprintf('  Min angle [deg] min/mean  : %.2f / %.2f\n', min(minAng),mean(minAng));
fprintf('  Frac. AR > 3             : %.1f%%\n', 100*mean(AR>3));
fprintf('  Frac. min-angle < 20 deg : %.1f%%\n', 100*mean(minAng<20));
fprintf('  Frac. min-angle < 10 deg : %.1f%%\n', 100*mean(minAng<10));
end
