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

%% Generate conformal surface sites using the UPR pipeline
% surfaceSites3D places equidistant site pairs on both sides of each
% triangulated surface so that the Voronoi faces align with the surfaces.
R = @(p) fGs * ones(size(p,1), 1);
F = surfaceSites3D({t1, t2}, {R, R});

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

% Remove reservoir sites that are too close to any surface site
% (sufficient conformity condition: no site may lie inside the balls
%  centred at the surface triangulation vertices)
rSites = surfaceSufCond3D(rSites, F.c.CC, F.c.R);

%% Merge surface sites and reservoir sites, then build grid
sites = [F.f.pts; rSites];

fprintf('Total number of sites = %d\n', size(sites,1));

G = mirroredPebi3D(sites, bdr);
G = computeGeometry(G);

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
end
