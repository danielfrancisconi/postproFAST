function EP(agir)
%=============================================================================================
%                     Entrega 05 - Exerc�cio Proposto Final                                     
%---------------------------------------------------------------------------------------------
% PNV5813 - M�todos Num�ricos Aplicados � Engenharia Naval e Oce�nica                        
% Ministrado pelo Prof. Dr. Andr� Lu�s Condino Fujarra                                       
%                                                                                            
%                                                             Universidade de S�o Paulo 2020 
%*********************************************************************************************
% Aluno:    Daniel Francisconi Oliveira                                                      
% NUSP:     8990622                                                                          
%*********************************************************************************************
%                                                           C�digo adaptado de (FUJARRA 2020)
%---------------------------------------------------------------------------------------------
% Fun��o para p�s tratamento de dados produzidos pelo programa FASTv8 para simula��es de 
% turbinas
%=============================================================================================

%.............................................................................................
% Par�metros/Vari�veis globais
global dados t titulo unidade filename legenda legenda2 sel grupo eixoX nomeX analise skip;
global Pai1 Pai2 Pai3 mensagem salva_auto Fig1 corrigido
%.............................................................................................
% Inst�ncia l�gica para prepara��o do ambiente de an�lise
if nargin==0
  clc, close all, clear all, more off, warning off
%   pkg load statistics % Necess�rio para o Boxplot
%   set(0,'defaultaxesfontname','helvetica')
%   set(0,'defaultaxesfontsize',20)
  corrigido = false; % corre��o aproamento j� feita -> true
  salva_auto = 0; % salvamento autom�tico
  EP('montamenu')
%.............................................................................................
% Inst�ncia l�gica para montagem do menubar
    elseif strcmp(agir,'montamenu')
      % Criando a janela
      Fig1      = figure('visible','off',...
          'name','P�s processamento de dados FASTv8',...
          'numbertitle','off',...
          'menubar','none',...
          'position',[100 200 500 400]);
      % Cria a 1� fam�lia da menubar
      Pai1        = uimenu(...
          'label','&Carregar',...
          'callback','EP(''carregar'')');
      % Cria a 2� fam�lia da menubar
      Pai2        = uimenu('enable','off',...
          'label','&Analisar',...
          'callback','EP(''def'')');
      % Cria a 3� fam�lia da menubar
      Pai3      = uimenu('label','&Reiniciar',...
          'accelerator','q',...
          'callback','EP(''reinicia'')');
      % Cria a 4� fam�lia da menubar
      Pai4      = uimenu('label','&Fechar');
      Filho1_Pai4 = uimenu(Pai4,...
          'label','&figuras',...
          'callback','EP(''limpa'')');
      Filho2_Pai4 = uimenu(Pai4,...
          'accelerator','q',...
          'label','tudo',...
          'callback','EP(''fechar'')');
      % Torna a figura vis�vel
      set(Fig1,'visible','on');
      title ('Ol�, carregue os dados para iniciar')
      texto = {'Esta � uma ferramenta para p�s processamento de dados',...
                'do FASTv8 para simula��es de turbinas e�licas. O', ...
                'programa foi desenvolvido como trabalho final do aluno'...
                'Daniel Francisconi Oliveira para a disciplina PNV5813'...
                'ministrada pelo Prof. Dr. Andr� Lu�s Condino Fujarra.'};
%                'Para utiliz�-lo basta que o output do FAST seja um .txt'}
      mensagem = uicontrol ('style','text','string',texto,...
      'units','normalized','position',[0.1 .2 .8 .2]);
      
      set(gca,'visible','off')
%.............................................................................................
% Inst�ncia l�gica para leitura do arquivo *.out output do FASTv8
elseif strcmp(agir,'carregar')
title('Fazendo a leitura dos arquivos ...')
% [Nome,Caminho,Fez] = uigetfile({'*.out', "Formato de Dados Esperado"},'MultiSelect', 'on');
  [filename,pathname,~] = uigetfile({'*.out', 'Formato de Dados Esperado'},'MultiSelect', 'on');
  % N�mero de linhas do cabe�alho:
  cab = 8;
  if isequal(pathname,0)
    % No file Selected
    disp('Nenhum arquivo selecionado!');
    return;
  elseif iscell(filename)
    % Multiple Files Selected
    FilesCount=length(filename);
    for loop1=1:FilesCount
      sFileName=filename{loop1};
      if ~isequal(sFileName,0)
        sFullFileName=fullfile(pathname,sFileName);
        fprintf ('leitura %d...', loop1)
        dados{loop1} = importdata(sFullFileName,'\t', cab);
        disp (' finalizada')
      end
    end
  else
    % One File Selected
    sFileName=filename;
    sFullFileName=fullfile(pathname,sFileName);
    %[~,sModelName,~] = fileparts(sModelName);
    fprintf('leitura...')
    dados{1} = importdata(sFullFileName, '\t', cab);
    disp(' finalizada')
  end
  % Cabe�alho (considerando que todas as simula��es tenham o mesmo):
  titulo = strsplit(char(dados{1}.textdata(cab-4)));% octave: cab-1
  titulo(length(titulo)) = []; % Deletar o �ltimo t�tulo (vazio)
  unidade = dados{1}.colheaders;
  % Vetores tempo:
  for loop2=1:length(dados)
    t{loop2}=dados{loop2}.data(1:end,1);
  end
  t = cell2mat(t); % tratar como cell ou matriz? Para usar mat, os vetores t devem ser iguais
  
  set(Pai1,'enable','off')
  set(Pai2,'enable','on')
  title('Clique em analisar e escolha quantos quiser')
  set(mensagem, 'visible', 'off')
  mensagem = uicontrol ('style','text','string','Selecione se deseja registrar os gr�ficos em .pdf',...
  'units','normalized','position',[0.1 .6 .8 .1]);
  mensagem = uicontrol ('style','checkbox','string','salvamento autom�tico',...
  'units','normalized','position',[0.3 .5 .4 .1], 'callback','EP(''auto'')');
%.............................................................................................
% Inst�ncia l�gica para definir salvamento autom�tico dos gr�ficos
elseif strcmp(agir,'auto')
if salva_auto
    salva_auto = false;
else
    salva_auto = true;
end
% falta cancelar o salvamento automatico
%.............................................................................................
% Inst�ncia l�gica para reiniciar
elseif strcmp(agir,'reinicia')
EP
%.............................................................................................
% Inst�ncia l�gica para definir plot
elseif strcmp(agir,'def')
  % Sele��o dos tipos de an�lise
  op_def1 = {'S�rie temporal', 'Transformada de Fourier', 'M�dia e Desvio', 'Boxplot', 'Aproamentos'}; 
  [analise,ok] = listdlg('Name','An�lises','PromptString','Que an�lise deseja fazer:',...
  'liststring',op_def1, 'SelectionMode', 'multiple');
  if ok==1
    disp('O usu�rio selecionou a(s) an�lise(s):')
    for u=1:numel(analise)
      disp(sprintf('\t%s', op_def1{analise(u)}))
    end
  else
    disp('O usu�rio desistiu da sele��o.')
    pause(1)
    % clc,clear all
    return
  end
% Sele��o de quais vari�veis plotar
if ~isempty(find(analise==1)) || ~isempty(find(analise==2))|| ~isempty(find(analise==3))|| ~isempty(find(analise==4)) % se todos exceto aproamentos
  op_def =titulo ;
  for loop4=1:length(titulo)
    op_def{loop4} = [num2str(loop4) '-' char(titulo(loop4))]; % enumerar op��es de titulo. octave: op_def(loop4).
  end
  [sel,ok] = listdlg('Name','Op��es','PromptString','Que vari�veis deseja plotar:',...
  'liststring',op_def, 'SelectionMode', 'multiple',...
  'InitialValue', [15 16 17 18 19 20]); % Op��es default
  if ok==1
    disp('O usu�rio selecionou as vari�veis:')
    for u=1:numel(sel)
      disp(sprintf('\t%s', op_def{sel(u)}))
    end
  else
    disp('O usu�rio desistiu da sele��o.')
    pause(1)
%    clc,clear all
    return
  end
end

if ~isempty(find(analise==1)) || ~isempty(find(analise==2))  % Se S�rie Temporal ou An�lise de Fourier
  % Defini��o da legenda
%   nome_leg = strsplit(num2str([1:1:length(dados)])); % Cria��o da legenda default (1,2,3...). octave:num2cell
%   nome_leg = {'0^o','10^o','20^o','30^o','40^o','50^o','60^o'}; % legenda defautlt para aproamentos
  nome_leg = {'7-0^o','7-10^o','7-20^o','7-30^o','7-40^o','7-50^o','7-60^o','12-0^o','12-10^o','12-20^o','12-30^o','12-40^o','12-50^o','12-60^o','20-0^o','20-10^o','20-20^o','20-30^o','20-40^o','20-50^o','20-60^o'};
  legenda = inputdlg(filename,'Defini��o da legenda',1, nome_leg); % Usu�rio pode escrever a legenda ou cancelar (sem legenda)
  if (isempty (legenda))
    legenda = 'off';
    disp ('O usu�rio desativou as legendas')
  else
    disp('O usu�rio definiu as legendas:')
    disp(legenda)
  end
end

% Defini��o do intervalo de an�lise
if ~isempty(find(analise==2)) || ~isempty(find(analise==3)) || ~isempty(find(analise==4))  % Se Fourier, M�dia/Desvio ou Boxplot 
  skip = {'2000'}; % Itera��o inicial para an�lise default
%   skip = {'150'}; % Voc� pode modificar o default aqui
  skip = inputdlg('Analisar a partir da itera��o n�mero:','Intervalo de An�lise',1, skip); % Usu�rio pode escrever a legenda ou cancelar (sem legenda)
  if (isempty (skip))
    disp('O usu�rio desistiu da sele��o.')
    pause(1)
    return
  else
    skip = str2num(cell2mat(skip));
    disp('O usu�rio definiu o intervalo de an�lise a partir da itera��o de n�mero:')
    disp(skip)
    if skip==0 
       disp('Aten��o! N�o h� itera��o de n�mero 0. A an�lise ser� feita a partir de 1.')
       skip =1;
    end
  end  
end

if ~isempty(find(analise==3)) || ~isempty(find(analise==4))  % Se M�dia/Desvio ou Boxplot
  %%% Agrupar os dados: (cada grupo deve ter o mesmo n�mero de dados)
  grupo = {}; % reiniciar para nova an�lise
  ngrupos= {'1'}; % N�mero de grupos default
  ngrupos = inputdlg('Defina o n�mero de grupos','N�mero de agrupamentos',1,ngrupos);
  ngrupos = str2num(cell2mat(ngrupos));
  qtd = idivide(length(dados),int32(ngrupos)); % Quantidade default de arquivos por grupo (divis�o inteira)
  for n=1:ngrupos
    grup_default = [((n-1)*qtd+1):1:n*qtd];
    [grupo{n},ok] = listdlg('Name','Definir Grupo','PromptString',[num2str(n) 'o grupo:'],...
    'liststring',filename, 'SelectionMode', 'multiple',...
    'InitialValue', grup_default); % Op��es default
    if ok==1
      disp(['Arquivos do ' num2str(n) 'o grupo:'])
      if iscell(filename) % se mais de um arquivo carregado
        for u=1:numel(grupo{n})
          disp(sprintf('\t%s', filename{grupo{n}(u)}))
        end
      else   % se apenas um arquivo carregado  
          disp(filename)
      end
    else
      disp('O usu�rio desistiu da sele��o.')
      pause(1)
    % clc,clear all
    return
    end
  end
  % Defini��o da legenda dos grupos
  nome_leg2 =  strsplit(num2str(1:1:ngrupos));%num2cell([1:1:ngrupos]) % Cria��o da legenda default (1,2,3...)
  for u=1:ngrupos
    nome_leg3{u} = ['Grupo ' num2str(u) ':'];
  end
  legenda2 = inputdlg(nome_leg3,'Defini��o da legenda',1,nome_leg2);%,1, nome_leg2); % Usu�rio pode escrever a legenda ou cancelar (sem legenda)
  if (isempty (legenda2))
    legenda2 = 'off';
    disp ('O usu�rio desativou as legendas')
  else
    disp('O usu�rio definiu as legendas:')
    disp(legenda2)
  end
  % Defini��o do t�tulo do eixo x:
%   nomeX = {'Velocidade de vento (m/s)'}; % Nome default
  nomeX = {'�ngulo de aproamento (deg)'};
  nomeX = inputdlg('T�tulo do eixo x:','Defini��o do eixo x',1, nomeX); % Usu�rio pode escrever a legenda ou cancelar (sem legenda)
  if (isempty (nomeX))
    disp ('O usu�rio desativou o t�tulo do eixo x')
  else
    disp('O usu�rio definiu o t�tulo do eixo x:')
    disp(nomeX)
  end
  % Defini��o das coordenadas do eixo x
  nome_leg2 =  strsplit(num2str(1:1:numel(grupo{1}))); % considerando que o num de pontos de cada grupo � igual do grupo{1}
  nome_leg2 = {'0' '10' '20' '30' '40' '50' '60'};
  eixoX = inputdlg(nome_leg2,'Coordenadas x',1,nome_leg2);%,1, nome_leg2); % Usu�rio pode escrever a legenda ou cancelar (sem legenda)
  if (isempty (eixoX))
    disp('O usu�rio desistiu da sele��o.')
    pause(1)
    return
  %  disp ('O usu�rio n�o definiu as coordenadas do eixo X para o plot das m�dias')
  %  disp ('O programa utilizar� a op��o default: [1,2,3..]')
  %  eixoX = nome_leg2;
  else
    disp('O usu�rio definiu as coordenadas do eixo X:')
    disp(eixoX)
  end
end %if
% Corre��o nos aproamentos
if ~isempty(find(analise==5))% Se corre��o de aproamentos
    Aproamentos();
end
EP('plot')
%-----------------
% Plotar
%------------------
elseif strcmp(agir,'plot')
%title('Clique em Reiniciar para uma nova an�lise')

  % Serie Temporal:
if ~isempty(find(analise==1))
  for u=1:numel(sel)
    coluna = sel(u);
    SerieTemporal(coluna)
  end
end

  % An�liise de Fourier:
if ~isempty(find(analise==2))
  for u=1:numel(sel)
    coluna = sel(u);
    EspectroAmplitude(coluna);
  end
end

  % M�dia e Desvio Padr�o:
if ~isempty(find(analise==3))
  for u=1:numel(sel)
    coluna = sel(u);
    MediaDesv(coluna)
  end
end
  
  % BoxPlot
if ~isempty(find(analise==4))
  for u=1:numel(sel)
    coluna = sel(u);
    BoxPlot(coluna)
  end
end
disp ('Pronto!')
%set(Pai2,'enable','off')
% set(mensagem, 'visible', 'off')
texto = {'Fim da an�lise',...
        'Voc� pode fazer uma nova an�lise com os mesmos arquivos',...
        'ou Reiniciar para carregar outros arquivos'};
mensagem = uicontrol (Fig1,'style','text','string',texto,...
  'units','normalized','position',[0.1 .6 .8 .2]);
%------------------------------------
% Inst�ncia l�gica para finaliza��o
%-------------------------------------
elseif strcmp(agir,'fechar')
 clc, close all, clear all
 disp('')
 disp('An�lise finalizada!')
%--------------------------------------
% Inst�ncia l�gica para fechar figuras
%--------------------------------------
% Fecha todas as figuras exceto o GUI
elseif strcmp(agir,'limpa')
 disp('Fechando as figuras...')
 set(Fig1,'HandleVisibility', 'off');
 close all
 set(Fig1,'HandleVisibility', 'on');
 disp('Pronto.')
end
end
%---------------------------------------------------
function SerieTemporal(coluna)
%---------------------------------------------------
%% Plota a s�rie temporal
global titulo unidade legenda dados t salva_auto
gtitulo = titulo{coluna};
gunidade = unidade(coluna);
disp(sprintf('plotando \t%s', gtitulo))
figure ('units', 'normalized', 'position', [.1 .08 .8 .8],'visible','on');
hold on
for loop3=1:length(dados)
  plot(t(:,loop3),dados{loop3}.data(1:end, coluna));
  title ([char(gtitulo), ' para diferentes �ngulos de aproamento  e U=7m/s'])
  xlabel('tempo (s)'), ylabel(gunidade), box on
%   set(gcf, 'visible','on')
  legend (legenda)
end
if salva_auto % se salvamento autom�tico ligado
  salvar(gcf, ['SerieTempo_' char(gtitulo)])
end
end
%---------------------------------------------------
function MediaDesv(coluna)
%---------------------------------------------------
%% Plota a m�dia e o desvio padr�o para cima e para baixo
global titulo unidade  dados grupo legenda2 eixoX nomeX skip salva_auto
% Prepar o sinal
sinal = [];
for loop6=1:length(dados) 
   sinal(:,loop6) = dados{loop6}.data(skip:end, coluna);
%    sinal{loop6} = dados{loop6}.data(skip:end, coluna);
end
  % C�lculo da m�dia:
  med = mean(sinal);
  % C�lculo do desvio padr�o:
  dev = std(sinal);
% Plot
gtitulo = titulo{coluna};
gunidade = unidade(coluna);
x = str2double(eixoX); %octave: cell2mat
fprintf('plotando M�dias: \t%s \n', gtitulo)
figure ('units', 'normalized', 'position', [.1 .08 .4 .5],'visible','on');
hold on
for n=1:numel(grupo) 
    medGrupo = med(:,grupo{n});
    devGrupo = dev(:,grupo{n});
    errorbar(x, medGrupo, devGrupo,'-s') % Plotar m�dia e desvio padr�o
%     plot(x, medGrupo,'-s') % Plotar somente a m�dia
end
hold off
title ([char(gtitulo), ' para diferentes aproamentos e velocidades de vento'])
xlabel(nomeX), ylabel(gunidade), box on, grid on
set(gca, 'fontsize',10)
legend(legenda2)
if salva_auto % se salvamento autom�tico ligado
  salvar(gcf, ['Media_' char(gtitulo)])
end
end
%---------------------------------------------------
function BoxPlot(coluna)
%---------------------------------------------------
%% Plota a mediana, os quartis e os outliers
global titulo unidade  dados grupo legenda2 eixoX nomeX skip salva_auto
sinal = [];
for u=1:length(dados) 
   sinal(:,u) = dados{u}.data(skip:end, coluna);
end

% Plot
gtitulo = titulo{coluna};
gunidade = unidade(coluna);
x = str2double(eixoX);
fprintf('plotando BoxPlot: \t%s \n', gtitulo)
figure ('units', 'normalized', 'position', [.1 .08 .5 .5],'visible','on'); % erro quando visible off: somente fig pares
hold on
for n=1:numel(grupo) 
  boxplot(sinal(:,grupo{n}), x,'OutlierSize',2,'Symbol', '.');
end
hold off
title ([char(gtitulo), ' para diferentes �ngulos de aproamento e velocidades de vento'])
xlabel(nomeX), ylabel(gunidade), box on
% grid on
set(gca, 'fontsize',10)
legend(legenda2)

if salva_auto % se salvamento autom�tico ligado
  salvar(gcf, ['BoxPlot_' char(gtitulo)])
end
end
%---------------------------------------------------
function [freq,Amp,fd,Ad]=EspectroAmplitude(coluna)
%---------------------------------------------------
%% Calcula o espectro de amplitude do sinal
%[freq,Amp,fd,Ad]=EspectroAmplitude(t,sinal)
% Entradas - t -> vetor de tempo
%            sinal
% saida    - freq -> vetor de frequencias [Hz]
%          - Amp  -> vetor de amplitudes
%          - fd   -> frequencia dominante
%          - Ad   -> amplitude na frequencia dominante
%% 
global dados legenda unidade t titulo skip salva_auto
% Prepar o sinal
sinal = [];
for loop6=1:length(dados) 
  sinal(:,loop6) = dados{loop6}.data(skip:end, coluna);
end
% FFT:
N=length(sinal); % antes length(t)
% deltat=t(2)-t(1);
% % deltaf=1/deltat;
% freq=[0:N-1]*deltaf/N;
Xs=fft(sinal);
Amp=abs(Xs)/N;
Amp = [Amp(1,:);2*Amp(2:N/2+1,:)];
freq = [0:N/2]'/t(N); % antes length(t)
%freq = repmat(freq,1,2);
%Amp=2*Amp(1:fix(N/2));
%freq=freq(1:fix(N/2));
[Ad,indice]=max(Amp(2:end,:)); % 2:end para evitar freq zero como maximo
fd=freq(indice);
Td = 1./fd;
% plot:
gtitulo = titulo{coluna};
gunidade = unidade(coluna);
figure ('units', 'normalized', 'position', [.1 .08 .5 .5],'visible','on');
fprintf ('Fr�quencias e per�odos dominantes em %s: \n ',gtitulo)
fprintf('%f Hz  %.1f s \n',[fd Td]')
loglog (freq,Amp)
title (['An�lise de Fourier de ', gtitulo])
legend (legenda)
xlabel('frequ�ncia (Hz)'), ylabel(gunidade)
set(gca, 'fontsize',10)
box on
grid on
if salva_auto % se salvamento autom�tico ligado
  salvar(gcf, ['Fourier_' char(gtitulo)])
end

end

%---------------------------------------------------
function Aproamentos()
%---------------------------------------------------
%% Modifica os sinais de surge sway roll e pitch fazendo a mudan�a de coordenadas 
% fazendo do novo eixo x a dire��o de incid�ncia de vento e onda (alinhados)
% os dados devem estar ordenados em velocidade 1 (aproamentos), velocidade
% 2 (aproamentos) ...
global dados corrigido
if corrigido
    disp('Corre��o no aproamento j� foi feita anteriormente')
else
disp('Fazendo a corre��o no aproamento')
colsurge = 15;
colsway = 16;
colroll = 18;
colpitch= 19;
for u=1:length(dados) 
    surge{u}=dados{u}.data(1:end,colsurge);
    sway{u}=dados{u}.data(1:end,colsway);
    roll{u}=dados{u}.data(1:end,colroll);
    pitch{u}=dados{u}.data(1:end,colpitch);
end
for n=1:1 % n�mero de velocidades
for u=1:7 % n�mero de aproamentos
    surge_novo{u+7*(n-1)} = surge{u+7*(n-1)}*cos(10*(u-1)*pi/180) + sway{u+7*(n-1)}*sin(10*(u-1)*pi/180);
    sway_novo{u+7*(n-1)} = sway{u+7*(n-1)}*cos(10*(u-1)*pi/180) - surge{u+7*(n-1)}*sin(10*(u-1)*pi/180);
    roll_novo{u+7*(n-1)} = roll{u+7*(n-1)}*cos(10*(u-1)*pi/180) + pitch{u+7*(n-1)}*sin(10*(u-1)*pi/180);
    pitch_novo{u+7*(n-1)} = pitch{u+7*(n-1)}*cos(10*(u-1)*pi/180) - roll{u+7*(n-1)}*sin(10*(u-1)*pi/180);
end
end
% Modificar sinal original:
for u=1:length(dados) 
    dados{u}.data(1:end,colsurge) = surge_novo{u};
    dados{u}.data(1:end,colsway) = sway_novo{u};
    dados{u}.data(1:end,colroll) = roll_novo{u};
    dados{u}.data(1:end,colpitch) = pitch_novo{u};
end %for
corrigido = true;
end %if
end %function
%---------------------------------------------------
function salvar(gcf, nome)
%---------------------------------------------------
%% Salvar uma figura
%%save to pdf
set(gcf,'Units','inches');
position = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0  position(3:4)],'PaperSize',[position(3:4)]);
print (gcf, nome, '-dpdf','-r0')
end
