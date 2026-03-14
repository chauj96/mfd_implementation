function G3 = buildMRSTMesh(Lx, Ly, H)
    
    fl = {[0.3 * Lx, 0.4 * Ly; 0.7 * Lx, 0.8 * Ly]};
    wl = {[0.3 * Lx, 0.4 * Ly; 0.7 * Lx, 0.8 * Ly]};
    
    G = compositePebiGrid2D([0.1, 0.1], [Lx, Ly], ...
        'faceConstraints', fl, ...
        'cellConstraints', wl, ...
        'FCFactor', 0.35, ...
        'CCFactor', 0.35, ...
        'mlqtMaxLevel', 2);
    
    G = computeGeometry(G);
    
    xy = G.nodes.coords;
    x = xy(:,1);
    y = xy(:,2);
    
    tol = 1e-10;
    isBnd = abs(x-0) < tol | abs(x-Lx) < tol | abs(y-0) < tol | abs(y-Ly) < tol;
    
    % band-based distortion
    dx = zeros(size(x));
    dy = zeros(size(y));
    
    % strip centers and half-widths
    yCenters = [2.0, 3.5, 5.0, 6.5, 8.0];
    bandHalfWidth = 0.5;   % thickness of each strip
    
    for k = 1:length(yCenters)
        yc = yCenters(k);
    
        inBand = abs(y - yc) <= bandHalfWidth;
    
        % smooth weight inside each strip
        s = zeros(size(y));
        s(inBand) = 0.5 * (1 + cos(pi*(y(inBand) - yc)/bandHalfWidth));
    
        % give each strip a slightly different motion
        if mod(k,3) == 1
            % mostly horizontal shift with slight waviness
            dx = dx + 1.4 * s .* (0.8 + 0.2 * sin(2*pi*x/Lx));
            dy = dy + 0.15 * s .* sin(pi*x/Lx);
    
        elseif mod(k,3) == 2
            % shear-like band
            dx = dx - 1.1 * s .* (x - 0.5*Lx)/Lx;
            dy = dy + 0.35 * s .* cos(2*pi*x/Lx);
    
        else
            % curved band
            dx = dx + 0.9 * s .* sin(1.5*pi*x/Lx);
            dy = dy - 0.25 * s;
        end
    end
    
    % boundary damping
    distBnd = min([x, Lx-x, y, Ly-y], [], 2);
    buffer = 1.2;
    mask = min(distBnd/buffer, 1).^2;
    
    dx = dx .* mask;
    dy = dy .* mask;
    
    % apply distortion
    amp2D = 0.07;
    xy_new = xy;
    xy_new(~isBnd,1) = xy(~isBnd,1) + amp2D*dx(~isBnd);
    xy_new(~isBnd,2) = xy(~isBnd,2) + amp2D*dy(~isBnd);
    
    G.nodes.coords = xy_new;
    G = computeGeometry(G);
    
    figure;
    plotGrid(G);
    axis equal tight;
    title('2D grid');
    
    nLayers = 14; % we can adjust number of layers (z direction)
    G3 = makeLayeredGrid(G, nLayers);
    G3 = computeGeometry(G3);
    
    zmin = min(G3.nodes.coords(:,3));
    zmax = max(G3.nodes.coords(:,3));
    G3.nodes.coords(:,3) = H * (G3.nodes.coords(:,3) - zmin) / (zmax - zmin);
    
    G3 = computeGeometry(G3);
    
    fprintf('3D cells = %d\n', G3.cells.num);
    
    figure;
    plotGrid(G3);
    view(3);
    axis equal tight;
    title('3D multilayer grid');
end

