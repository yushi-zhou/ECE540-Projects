clear all;
close all;
clc;

a = 0.033; % WG width
b = 0.014; % WG height
hmax = a/200; % Mesh size

% Use PDE Toolbox to create mesh
model = createpde();
gd = [3;        % Rectangle
      4;        % 4 points
      0; a; a; 0;  % x coordinates
      0; 0; b; b]; % y coordinates
sf = 'R1'; % Shape formula
ns = char('R1')';
dl = decsg(gd, sf, ns);
geometryFromEdges(model, dl);
generateMesh(model, 'Hmax', hmax, 'GeometricOrder', 'linear');

p = model.Mesh.Nodes;   % Node coordinates (2 x numNodes)
t = model.Mesh.Elements;     % Element connectivity (3 x numElements)
numNodes = size(p, 2);    
numElements = size(t, 2);   

figure;
pdemesh(model);
axis equal;
title('Triangular Mesh');
xlabel('x (m)');
ylabel('y (m)');

%Matrix Assembly
K = sparse(numNodes, numNodes); 
M = sparse(numNodes, numNodes); 

for e = 1:numElements
    nodes = t(:, e);         % 3x1 node indices
    p_e = p(:, nodes);       % 2x3 node coordinates
    x = p_e(1, :);           
    y = p_e(2, :); 
    
    A_e = abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))) / 2;
    
    b_i = [y(2) - y(3); x(3) - x(2)]; % Gradient coefficients
    b_j = [y(3) - y(1); x(1) - x(3)];
    b_k = [y(1) - y(2); x(2) - x(1)];
    
    % Local K
    K_e = (1 / (4 * A_e)) * [
        dot(b_i, b_i), dot(b_i, b_j), dot(b_i, b_k);
        dot(b_j, b_i), dot(b_j, b_j), dot(b_j, b_k);
        dot(b_k, b_i), dot(b_k, b_j), dot(b_k, b_k)
    ];
    
    % Local M
    M_e = (A_e / 12) * [
        2, 1, 1;
        1, 2, 1;
        1, 1, 2
    ];
    
    % Golbal
    for i = 1:3
        for j = 1:3
            K(nodes(i), nodes(j)) = K(nodes(i), nodes(j)) + K_e(i, j);
            M(nodes(i), nodes(j)) = M(nodes(i), nodes(j)) + M_e(i, j);
        end
    end
end

% Solve
numModes = 6; 
[V_all, D_all] = eigs(K, M, numModes+1, 'smallestabs'); 
V = V_all(:, 2:numModes+1);
D = D_all(2:numModes+1, 2:numModes+1);
k_c_squared = diag(D); 
k_c = sqrt(k_c_squared); 
c = 3e8; 
cutoff_freq = k_c * c / (2 * pi); 

fprintf('Rectangular waveguide dimensions: %.2f mm x %.2f mm\n', a*1000, b*1000);
fprintf('\nTE Mode Analysis Results:\n');
fprintf('-------------------------------------------\n');
fprintf('Mode  |  Cutoff k(rad/m)  |  Cutoff f(GHz)\n');
fprintf('-------------------------------------------\n');
for i = 1:numModes
    fprintf(' %d    |     %8.2f      |     %8.2f\n', ...
        i, k_c(i), cutoff_freq(i)/1e9);
end


%Theoretical
m_values = [1, 2, 0, 1, 3, 2]; % For TE10, TE01, TE11, TE20, TE02, etc.
n_values = [0, 0, 1, 1, 0, 1];
theoretical_kc = zeros(size(m_values));
theoretical_fc = zeros(size(m_values));
for i = 1:length(m_values)
    m = m_values(i);
    n = n_values(i);
    theoretical_kc(i) = sqrt((m * pi / a)^2 + (n * pi / b)^2);
    theoretical_fc(i) = theoretical_kc(i) * c / (2 * pi);
end

fprintf('\nComparison with Theoretical Values:\n');
fprintf('-------------------------------------------\n');
fprintf('Mode   |  Theory f_c(GHz)  |  FEM f_c(GHz)  |  Error(%%)\n');
fprintf('-------------------------------------------\n');

total_error = 0;
for i = 1:numModes
    
    % Find closest theoretical mode
    [~, idx] = min(abs(theoretical_fc - cutoff_freq(i)));
    error_percent = abs(cutoff_freq(i) - theoretical_fc(idx)) / theoretical_fc(idx) * 100;
    m = m_values(idx);
    fprintf('TE_%d%d  |     %8.2f    |     %8.2f     |   %6.4f\n', ...
        m, n_values(idx), theoretical_fc(idx)/1e9, cutoff_freq(i)/1e9, error_percent);
    total_error = total_error + error_percent;
end
fprintf('average error: %.4f%%\n', total_error / numModes);


%Plot
figure;
nRows = ceil(sqrt(numModes));
nCols = ceil(numModes / nRows);
mode_names = cell(1, numModes);

for i = 1:numModes
    [~, idx] = min(abs(theoretical_fc - cutoff_freq(i)));
    m = m_values(idx);
    n = n_values(idx);
    mode_names{i} = sprintf('TE_{%d%d}', m, n);
end

for i = 1:numModes
    subplot(nRows, nCols, i);
    field = -V(:, i);
    trisurf(t(1:3,:)', p(1,:), p(2,:), field, 'EdgeColor', 'none');
    view(2); % 2D view
    axis equal tight;
    title(sprintf('%s - FEM field f_c = %.2f GHz', mode_names{i}, cutoff_freq(i)/1e9));
    xlabel('x (m)');
    ylabel('y (m)');
    colorbar;
    colormap(jet);
end
sgtitle('Metallic Waveguide TE Modes - FEM Solution');

%plot Theoretical
figure;
[X, Y] = meshgrid(linspace(0, a, 100), linspace(0, b, 100));
%mode_names = cell(1, length(m_values));
mode_names = {'TE_{10}', 'TE_{20}', 'TE_{01}', 'TE_{11}', 'TE_{30}', 'TE_{21}'};

for i = 1:length(m_values)
    m = m_values(i);
    n = n_values(i);

    if n == 0
        H_z = cos(m * pi * X / a);
    elseif m == 0
        H_z = cos(n * pi * Y / b);
    else
        H_z = cos(m * pi * X / a) .* cos(n * pi * Y / b);
    end

    subplot(nRows, nCols, i);
    surf(X, Y, H_z, 'EdgeColor', 'none');
    view(2); % 2D view
    axis equal tight;
    title(sprintf('%s - Theory f_c = %.2f GHz', mode_names{i}, theoretical_fc(i)/1e9));
    xlabel('x (m)');
    ylabel('y (m)');
    colorbar;
    colormap(jet);
end
sgtitle('Rectangular Waveguide TE Modes - Theoretical Solution');