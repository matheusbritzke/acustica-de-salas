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
            fn_m{modo} = (sprintf('%.f, %.f, %.f', nx, ny, nz));
            modo = modo+1;
       end
    end
end


fn = fn(2:end);
fn_m = fn_m(2:end);

figure;
stem(fn, ones(length(fn)), 'Color', 'r', 'LineWidth', 1.2, 'Marker', 'none')
ylim([0,1.4])
yticks(1)
xlim([f(1),f(end)])

xline(fn(1),'--k','LineWidth',1.3)
xline(Fs,'--k','LineWidth',1.3)
title('Modos da sala retangular','FontSize', 14);
subtitle (sprintf('Dimensões: %.2fm x %.2fm x %.2fm',Lx,Ly,Lz), 'FontSize', 12)
ylabel('Núemro de modos [-]', 'FontSize', 12)
xlabel('Frequência [Hz]', 'FontSize', 12)
pbaspect([10 6 1]);
grid on

%ajustar eixo x
xticks([20,40,60,80,100,200,300,400])
ax = gca; % Obtém o handle do eixo atual
ax.XScale = 'log'; % Define a escala X como logarítmica

%regiões

text(f(1) , 1.3, ' Região X','HorizontalAlignment', 'left', 'FontSize',12)
text(fn(1) , 1.3, ' Região A','HorizontalAlignment', 'left', 'FontSize',12)
text( Fs, 1.3, ' Região B','HorizontalAlignment', 'left', 'FontSize',12)

% Ajustar propriedades para exportação
set(gcf, 'PaperUnits', 'centimeters');      
set(gcf, 'PaperSize', [45, 45*6/10]);            % Tamanho da página do PDF
set(gcf, 'PaperPosition', [0, 0, 45, 45*6/10]);   % Ocupa toda a página

% num = 1:1:length(fn);
% 
% 
% T= table (fn_m(:), fn(:), 'VariableNames',{'Modos','f_n'});
% T = sortrows (T, 'f_n');
% num = table (num(:), 'VariableNames',{'nº'});
% T = [num,T];
% disp (T)
% 
% writetable(T, 'teste', ...
%                'Delimiter', ';', ...
%                'WriteRowNames', true, ...
%                'WriteVariableNames', true);