% =========================================================================
% Script: formas_modais.m
%
% Descrição:
%   Obejtivo da etapa 3 do projeto final de acústica de salas
%
%
%
% Autor: Matheus Britzke
% Data: 19/05/2025
% =========================================================================

close all; clear; clc;



%parametros do ar
c0 = 343; %[m/s] velocidade do som

% parâmetros da sala:

Lx= 4.37; % [m] maior dimensão da sala
Ly= 3.36; % [m] segunda dimensão no horizontal
Lz= 2.85; % [m] altura

Area = Lx*Ly; % area
V = Area*Lz; % [m^3] volume da sala

T60 = 0.5; % [s] 0.5 pq é o pior cenário
Fs = 2000*sqrt(0.5/V); % [Hz] frequ~encia de shroeder

f=20:1:round(Fs)+200; % [Hz] espectro de frequências



% FORMA MODAL
% qual modo calcular?
nx = 3;
ny = 0;
nz = 2;

fM = (c0/2)*sqrt( (nx/Lx)^2 + (ny/Ly)^2 + (nz/Lz)^2 );

% formando os eixos
N = Lx*100; %resolução
eixo_x = linspace(0, Lx, N);
eixo_y = linspace(0, Ly, N);
eixo_z = linspace(0, Lz, N);

modo_x = nx > 0;
modo_y = ny > 0;
modo_z = nz > 0;

dim_ativa = modo_x + modo_y + modo_z;

set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12)
switch dim_ativa
    case 1  % Modo 1D
        if modo_x
            [X, Y] = meshgrid(eixo_x, eixo_y);
            phi = cos(nx*pi*X/Lx);  % Repete valor de phi em Y
            imagesc(eixo_x, eixo_y, abs(phi));
            xlabel('x (m)', 'FontSize',12); ylabel('y (m)', 'FontSize',12);
            c = colorbar;
            c.Label.String = '|\psi(x)|';
        elseif modo_y
            [X, Y] = meshgrid(eixo_x, eixo_y);
            phi = cos(ny*pi*Y/Ly);  % Repete valor de phi em X
            imagesc(eixo_x, eixo_y, abs(phi));
            xlabel('x (m)', 'FontSize',12); ylabel('y (m)', 'FontSize',12);
            c = colorbar;
            c.Label.String = '|\psi(y)|';
        elseif modo_z
            [X, Z] = meshgrid(eixo_x, eixo_z);
            phi = cos(nz*pi*Z/Lz);  % Repete valor de phi em X
            imagesc(eixo_x, eixo_z, abs(phi));
            xlabel('x (m)', 'FontSize',12); ylabel('z (m)','FontSize',12);
             c = colorbar;
             c.FontSize = 13;
             c.Label.String = '|\psi(z)|';
        end

        colormap(gray);
        axis xy;
        title(sprintf('Modo axial [%.f, %.f, %.f]',nx,ny,nz),"FontSize",14)
        subtitle(sprintf('Frequência: %.2f Hz', fM),"FontSize",13)
        c.Ticks = [0.001,0.2,0.4,0.6,0.8,1];
        c.TickLabels = [0,0.2,0.4,0.6,0.8,1];
        c.Label.FontSize = 14;

    case 2  % Modo 2D
        if ~modo_z  % xy
            [X, Y] = meshgrid(eixo_x, eixo_y);
            phi = cos(nx*pi*X/Lx) .* cos(ny*pi*Y/Ly);
            imagesc(eixo_x, eixo_y, abs(phi));
            xlabel('x (m)', 'FontSize',12); ylabel('y (m)', 'FontSize',12);
            c = colorbar;
            c.Label.String = '|\psi(x,y)|';
        elseif ~modo_y  % xz
            [X, Z] = meshgrid(eixo_x, eixo_z);
            phi = cos(nx*pi*X/Lx) .* cos(nz*pi*Z/Lz);
            imagesc(eixo_x, eixo_z, abs(phi));
            xlabel('x (m)', 'FontSize',12); ylabel('z (m)', 'FontSize',12);
            c = colorbar;
            c.Label.String = '|\psi(x,z)|';
        else  % yz
            [Y, Z] = meshgrid(eixo_y, eixo_z);
            phi = cos(ny*pi*Y/Ly) .* cos(nz*pi*Z/Lz);
            imagesc(eixo_y, eixo_z, abs(phi));
            xlabel('y (m)', 'FontSize',12); ylabel('z (m)', 'FontSize',12);
            c = colorbar;
            c.Label.String = '|\psi(y,z)|';
            
        end
        colormap(gray); axis xy;
        c.Label.FontSize = 14;
        c.Ticks = [0.001,0.2,0.4,0.6,0.8,1];
        c.TickLabels = [0,0.2,0.4,0.6,0.8,1];
        title(sprintf('Modo tangencial [%.f, %.f, %.f]',nx,ny,nz),"FontSize",14)
        subtitle(sprintf('Frequência: %.2f Hz', fM),"FontSize",13)
        hold off

    case 3  % Modo 3D
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
        figure;
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


        colormap(gray);
        c = colorbar;
        c.Label.FontSize = 14;
        c.Label.String = '|\psi(x,y,z)|';
        c.Ticks = [0.001,0.2,0.4,0.6,0.8,1];
        c.TickLabels = [0,0.2,0.4,0.6,0.8,1];
        xlabel('x (m)', 'FontSize',12); ylabel('y (m)', 'FontSize',12); zlabel('z (m)', 'FontSize',12);
        title(sprintf('Modo obllíquo [%.f, %.f, %.f]',nx,ny,nz),"FontSize",14)
        subtitle(sprintf('Frequência: %.2f Hz', fM),"FontSize",13)
        % view(30, 25);
        view(-37.5, 30);
        axis([0, Lx, 0, Ly, 0, Lz])
        grid on
        hold off
end

 pbaspect([Lx Ly Lz]);
 set(gcf, 'PaperUnits', 'centimeters');      
 set(gcf, 'PaperSize', [45, 45]);            % Tamanho da página do PDF
 set(gcf, 'PaperPosition', [0, 0, 45, 45]);   % Ocupa toda a página
