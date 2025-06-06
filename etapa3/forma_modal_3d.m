close all; clear; clc;

% Parâmetros do ar
c0 = 343; %[m/s] velocidade do som

% Parâmetros da sala:
Lx = 4.37; % [m] maior dimensão da sala
Ly = 3.36; % [m] segunda dimensão no horizontal
Lz = 2.85; % [m] altura

Area = Lx*Ly; % area
V = Area*Lz; % [m^3] volume da sala

T60 = 0.5; % [s] 0.5 pq é o pior cenário
Fs_shroeder = 2000*sqrt(0.5/V); % [Hz] frequência de Schroeder (renomeada para clareza)

% f = 20:1:round(Fs_shroeder)+200; % [Hz] espectro de frequências (não usado neste plot específico)

% FORMA MODAL
% qual modo calcular? (Altere estes valores para ver diferentes modos)
nx = 1;
ny = 1;
nz = 1;

fM = (c0/2)*sqrt( (nx/Lx)^2 + (ny/Ly)^2 + (nz/Lz)^2 );

% Formando os eixos
N_resolucao = 100; %resolução (renomeado para clareza)
eixo_x = linspace(0, Lx, N_resolucao);
eixo_y = linspace(0, Ly, N_resolucao);
eixo_z = linspace(0, Lz, N_resolucao);

% --- Geração das malhas e cálculo de phi para cada face ---

% Face xy no z=0 (Piso)
[MESH_X_xy0, MESH_Y_xy0] = meshgrid(eixo_x, eixo_y);
MESH_Z_xy0 = zeros(size(MESH_X_xy0));
phi_xy0 = cos(nx*pi*MESH_X_xy0/Lx) .* cos(ny*pi*MESH_Y_xy0/Ly) .* cos(nz*pi*MESH_Z_xy0/Lz);

% Face xy no z=Lz (Teto)
[MESH_X_xyLz, MESH_Y_xyLz] = meshgrid(eixo_x, eixo_y);
MESH_Z_xyLz = Lz*ones(size(MESH_X_xyLz));
phi_xyLz = cos(nx*pi*MESH_X_xyLz/Lx) .* cos(ny*pi*MESH_Y_xyLz/Ly) .* cos(nz*pi*MESH_Z_xyLz/Lz);

% Face xz no y=0 (Parede Lateral 1)
[MESH_X_xz0, MESH_Z_xz0] = meshgrid(eixo_x, eixo_z);
MESH_Y_xz0 = zeros(size(MESH_X_xz0));
phi_xz0 = cos(nx*pi*MESH_X_xz0/Lx) .* cos(ny*pi*MESH_Y_xz0/Ly) .* cos(nz*pi*MESH_Z_xz0/Lz);

% Face xz no y=Ly (Parede Lateral 2)
[MESH_X_xzLy, MESH_Z_xzLy] = meshgrid(eixo_x, eixo_z);
MESH_Y_xzLy = Ly*ones(size(MESH_X_xzLy));
phi_xzLy = cos(nx*pi*MESH_X_xzLy/Lx) .* cos(ny*pi*MESH_Y_xzLy/Ly) .* cos(nz*pi*MESH_Z_xzLy/Lz);

% Face yz no x=0 (Parede Frontal)
[MESH_Y_yz0, MESH_Z_yz0] = meshgrid(eixo_y, eixo_z); % Note a ordem para meshgrid
MESH_X_yz0 = zeros(size(MESH_Y_yz0));
phi_yz0 = cos(nx*pi*MESH_X_yz0/Lx) .* cos(ny*pi*MESH_Y_yz0/Ly) .* cos(nz*pi*MESH_Z_yz0/Lz);

% Face yz no x=Lx (Parede Traseira)
[MESH_Y_yzLx, MESH_Z_yzLx] = meshgrid(eixo_y, eixo_z); % Note a ordem para meshgrid
MESH_X_yzLx = Lx*ones(size(MESH_Y_yzLx));
phi_yzLx = cos(nx*pi*MESH_X_yzLx/Lx) .* cos(ny*pi*MESH_Y_yzLx/Ly) .* cos(nz*pi*MESH_Z_yzLx/Lz);

% --- Plot das superfícies 3D ---
fig = figure;
hold on;

% Plotando com a ordem correta (X, Y, Z) para o surf
surf(MESH_X_xy0,  MESH_Y_xy0,  MESH_Z_xy0,  abs(phi_xy0),  'EdgeColor', 'none');
surf(MESH_X_xyLz, MESH_Y_xyLz, MESH_Z_xyLz, abs(phi_xyLz), 'EdgeColor', 'none');

surf(MESH_X_xz0,  MESH_Y_xz0,  MESH_Z_xz0,  abs(phi_xz0),  'EdgeColor', 'none');
surf(MESH_X_xzLy, MESH_Y_xzLy, MESH_Z_xzLy, abs(phi_xzLy), 'EdgeColor', 'none');

surf(MESH_X_yz0,  MESH_Y_yz0,  MESH_Z_yz0,  abs(phi_yz0),  'EdgeColor', 'none');
surf(MESH_X_yzLx, MESH_Y_yzLx, MESH_Z_yzLx, abs(phi_yzLx), 'EdgeColor', 'none');



% Arestas da base (z=0)
plot3([0, Lx], [0, 0],   [0, 0],   'k', 'LineWidth', .5); % V1-V2
plot3([Lx, Lx], [0, Ly],  [0, 0],   'k', 'LineWidth', .5); % V2-V3
plot3([Lx, 0],  [Ly, Ly], [0, 0],   'k', 'LineWidth', .5); % V3-V4
plot3([0, 0],   [Ly, 0],  [0, 0],   'k', 'LineWidth', .5); % V4-V1

% Arestas do topo (z=Lz)
plot3([0, Lx], [0, 0],   [Lz, Lz], 'k', 'LineWidth', .5); % V5-V6
plot3([Lx, Lx], [0, Ly],  [Lz, Lz], 'k', 'LineWidth', .5); % V6-V7
plot3([Lx, 0],  [Ly, Ly], [Lz, Lz], 'k', 'LineWidth', .5); % V7-V8
plot3([0, 0],   [Ly, 0],  [Lz, Lz], 'k', 'LineWidth', .5); % V8-V5

% Arestas verticais
plot3([0, 0],   [0, 0],   [0, Lz],  'k', 'LineWidth', .5); % V1-V5
plot3([Lx, Lx], [0, 0],   [0, Lz],  'k', 'LineWidth', .5); % V2-V6
plot3([Lx, Lx], [Ly, Ly], [0, Lz],  'k', 'LineWidth', .5); % V3-V7
plot3([0, 0],   [Ly, Ly], [0, Lz],  'k', 'LineWidth', .5); % V4-V8

c = colorbar;
c.Ticks = [0.001,0.2,0.4,0.6,0.8,1];
c.TickLabels = [0,0.2,0.4,0.6,0.8,1];
c.Label.String = '|\psi(x,y,z)|';
c.Label.FontSize = 14;

hold off;

% --- Configurações do Gráfico ---
colormap(gray); % Ou 'jet', 'parula', etc.



% Rótulos dos eixos
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

% Título e Subtítulo
str_tipo_modo = '';
if nx~=0 && ny==0 && nz==0 || nx==0 && ny~=0 && nz==0 || nx==0 && ny==0 && nz~=0
    str_tipo_modo = 'Axial';
elseif (nx~=0 && ny~=0 && nz==0) || (nx~=0 && ny==0 && nz~=0) || (nx==0 && ny~=0 && nz~=0)
    str_tipo_modo = 'Tangencial';
elseif nx~=0 && ny~=0 && nz~=0
    str_tipo_modo = 'Oblíquo';
else % (0,0,0) - Não é um modo ressonante, mas pode ser plotado
    str_tipo_modo = 'Nulo (0,0,0)';
end
title(sprintf('Modo %s [%.f, %.f, %.f]', str_tipo_modo, nx, ny, nz), "FontSize", 14);
subtitle(sprintf('Frequência: %.2f Hz', fM), "FontSize", 13);

% Visualização e Eixos
view(-37.5, 30); % Ajuste o ângulo de visualização (azimute, elevação)
% lighting gouraud;
% camlight headlight; % Adiciona uma luz vinda da câmera
axis equal; % Garante que as proporções da sala sejam visualmente corretas
grid on;
box on; % Desenha uma caixa ao redor dos eixos

ax = gca;
ax.FontSize = 12; % Tamanho da fonte dos ticks dos eixos


% Ou mais simples:
xlim([0, Lx]); ylim([0, Ly]); zlim([0, Lz]); % Se não quiser margens
xticks([0,1,2,3,4,Lx]); xticklabels({'0','1','2','3','4', 'Lx'})
yticks([0,1,2,3,Ly]); yticklabels({'0','1','2','3', 'Ly'});
zticks([0,1,2,Lz]); zticklabels({'0','1','2','Lz'});

% Configurações para salvar em PDF (opcional)
pbaspect([Lx Ly Lz]); % Define a proporção dos eixos de dados para corresponder à sala
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20, 15]); % Ajuste o tamanho do papel
set(gcf, 'PaperPosition', [0, 0, 20, 15]); % Ajuste a posição no papel