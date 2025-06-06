close all; clear; clc;

c0 = 343; %[m/s] velocidade do som

% CERTO
Lx= 4.37; % [m] maior dimensão da sala
Ly= 3.36; % [m] segunda dimensão no horizontal
Lz= 2.85; % [m] altura

Area = Lx*Ly; % area
V = Area*Lz; % [m^3] volume da sala

T60 = 0.5; % [s] 0.5 pq é o pior cenário
Fs = 2000*sqrt(0.5/V); % [Hz] frequ~encia de shroeder

f=20:1:round(Fs)*2; % [Hz] espectro de frequências

ndm = 15 ; % numero final em cada eixo dos modos medidos
modo = 1; % um contador de indice
for nz=0:ndm
    for ny=0:ndm
       for nx=0:ndm
            fn(modo) = (c0/2)*sqrt( (nx/Lx)^2 + (ny/Ly)^2 + (nz/Lz)^2 );
            modo = modo+1;
       end
    end
end


fn = fn(2:end);
fn = sort(fn);
idx_valid_fn = fn >= f(1) & fn <= f(end);
fn = fn(idx_valid_fn);


function [Zc, kc]=DB(f,rho0,c0,sigma)
% Impedância característica do material
Zc = (rho0 * c0).* ( (1 + 9.08 .* ((1e3.*f)./sigma).^(-0.75)) - 1j * 11.9 .* ((1e3.*f)./sigma).^(-0.73) );

k0 = (2*pi.*f)./c0;
% Número de onda característico:
kc = k0.* ( (1 + 10.8.*((1e3.*f)./sigma).^(-0.70)) - 1j * (10.3 .* ((1e3.*f)./sigma).^(-0.59)) );
end

function alpha=alpha_membrana(f,sigma,D,d,dm,rhom)
omega = 2*pi*f; % [rad/s] Vetor de frequencia angular

% PROPRIEDADES AR
rho0 = 1.21; % [kg/m^3] densidade do ar
c0 = 343; % [m/s] velocidade do som no ar
k0 = omega/c0; % [rad/m] numero de onda

denssup = dm.*rhom; % [kg/m^2] Densidade superficial da membrana

Zm  = 1j .* omega .* denssup; % Impedância da membrana
[Zp, kp] = DB(f, rho0, c0, sigma); % Impedância característica do material poroso
Zsp = Zp .* ((cosh(kp.*d))./(sinh(kp.*d))); % Impedância no topo da camada de material poroso         
Zsi1 = (-1j.*Zsp.*rho0.*c0.*(1./tanh(k0.*(D-d)))) + ((rho0.*c0).^2);  
Zsi2 = Zsp - (1j.*rho0.*c0.*(1./tanh(k0.*(D-d))));

Zsi= Zsi1 ./ Zsi2 ; % Impedância no topo da camada de ar

Zs = Zm + Zsi; % Impedância de superfície do absorvedor

% Coeficiente de absorção alpha                   
reflex = (Zs - (rho0*c0)) ./ (Zs + (rho0*c0)); % Coeficiente de reflexão
alpha = 1 - (abs(reflex).^2);
end

% sigma = 25000; % [Ns/m^4] Resistividade ao fluxo do material poroso
% D = .20; % [m] comprimento total da cavidade do absorvedor
% d = 4e-2; % [m] Espessura do material poroso
% dm   = 25e-4; % [m] Espessura da membrana
% rhom = 5e3; % [kg/m^2] Densidade do material da membrana

% ordem dos dados: sigma, D, d, dm, rhom
mem= [
    10000, .40, 5e-2, 24e-4, 2.7e3;
    10000, .20, 5e-2, 10e-4, 1e3;
    ];



figure;
hold on
stem(fn, ones(length(fn)),'y', 'LineWidth', 1.2, 'Marker', 'none', 'HandleVisibility', 'off')
ylim([0,1.4])
yticks([0,0.2,0.4,0.6,0.8,1])
xlim([f(1),f(end)])

for i = 1:size(mem,1)
    cores_plot = lines(size(mem,1)); 
    sigma_atual = mem(i, 1);
    D_atual     = mem(i, 2);
    d_atual     = mem(i, 3);
    dm_atual    = mem(i, 4);
    rhom_atual  = mem(i, 5);
    nome_legenda = sprintf('Membrana %.f',i);
    alpha = alpha_membrana(f, sigma_atual, D_atual, d_atual, dm_atual, rhom_atual);
    semilogx(f, alpha, 'LineWidth', 3, 'DisplayName', nome_legenda, 'Color', cores_plot(i,:));
    
end
% Configurações finais do gráfico
ax=gca;
set(ax, 'XScale', 'log')
set(ax, 'TickLabelInterpreter', 'tex') 
xlabel('Frequência (Hz)')
legend('FontSize', 12)
xline(fn(1),'--k','LineWidth',1.3, 'HandleVisibility', 'off')
xline(Fs,'--k','LineWidth',1.3, 'HandleVisibility', 'off')
title('Modos e Absorção','FontSize', 14);
subtitle (sprintf('Dimensões: %.2fm x %.2fm x %.2fm',Lx,Ly,Lz), 'FontSize', 12)
ylabel('\alpha (Coeficiente de Absorção)', 'FontSize', 12)
xlabel('Frequência [Hz]', 'FontSize', 12)
pbaspect([10 8 1]);
grid on

hold off
%ajustar eixo x
xticks([20,40,60,80,100,200,300,400])
ax.XScale = 'log'; % Define a escala X como logarítmica
ax.FontSize = 12;
%regiões

text(f(1) , 1.3, ' Região X','HorizontalAlignment', 'left', 'FontSize',13)
text(fn(1) , 1.3, ' Região A','HorizontalAlignment', 'left', 'FontSize',13)
text( Fs, 1.3, ' Região B','HorizontalAlignment', 'left', 'FontSize',13)


% Ajustar propriedades para exportação
set(gcf, 'PaperUnits', 'centimeters');      
set(gcf, 'PaperSize', [45, 45*8/10]);            % Tamanho da página do PDF
set(gcf, 'PaperPosition', [0, 0, 45, 45*8/10]);   % Ocupa toda a página