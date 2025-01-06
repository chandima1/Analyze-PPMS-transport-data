clc
% clear all
close all

% Define the file path
filename1 = 'E:\Transport\CEIQW240911_B and CEQW240715_B\RvsT300K_1_8K_InSbBi240911_B_Chan3_InSb240715_B_Chan2.dat.csv'; %%InSb QW without cap - 240715-B
filename2 = 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\RvsT300K_1_8K_Chan2_CEQW241010_Chan3_CEQW241009.csv'; %%InSbBi QW without cap - 241009
filename3 = 'E:\Transport\CEIQW240911_B and CEQW240715_B\Symmetrized_Rxx_Data_CEQW240715_B.csv';
filename4 = 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\Symmetrized_Rxx_9T_Chan2_CEQW241010_Chan3_CEQW241009.csv';
filename5 = 'E:\Transport\CEIQW240911_B and CEQW240715_B\Symmetrized_Rxy_9T.csv';
filename6 = 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\Symmetrized_Rxy_9T_Chan2_CEQW241010_Chan3_CEQW241009.csv';





Data1 = readtable(filename1);
Data2 = readtable(filename2);
Data3 = readtable(filename3);
Data4 = readtable(filename4);
Data5 = readtable(filename5);
Data6 = readtable(filename6);

%% R vs T
figure1 = figure('WindowState', 'maximized');
plot(Data1{:,1},Data1{:,4},'b.-','MarkerSize',14,'LineWidth',2);
axis square
xlabel('Temperature (K)')
ylabel('Resistance (Ω)') ;
set(gca,'FontSize',22,'FontName','sans serif')
ax = gca;
xlim([0 302])
ylim([0 max(Data1{:,4})+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
set(gca,'color','white')
% exportgraphics(gca,'R vs T_InSb240715_B.tiff','Resolution',600)

figure2 = figure('WindowState', 'maximized');
plot(Data2{:,1},Data2{:,5},'b.-','MarkerSize',14,'LineWidth',2);
axis square
xlabel('Temperature (K)')
ylabel('Resistance (Ω)') ;
set(gca,'FontSize',22,'FontName','sans serif')
ax = gca;
xlim([0 302])
ylim([0 max(Data2{:,5})+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
set(gca,'color','white')
% exportgraphics(gca,'R vs T_CEQW241009.tiff','Resolution',600)

figure3 = figure('WindowState', 'maximized');
p1=plot(Data1{:,1},Data1{:,4},'r.-','MarkerSize',14,'LineWidth',2);
hold on
p2=plot(Data2{:,1},Data2{:,5},'b.-','MarkerSize',14,'LineWidth',2);
hold off
axis square
xlabel('Temperature (K)')
ylabel('Resistance (Ω)') ;
set(gca,'FontSize',22,'FontName','sans serif')
ax = gca;
xlim([0 302])
ylim([0 max(Data1{:,4})+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
lgd1=legend([p1,p2],'InSb','InSbBi','EdgeColor','none','FontSize',22,'FontName','sans serif', 'Position',[0.591623262162839 0.751885203413262 0.107682293442388 0.151152193185668]);
lgd1.Color = 'none';
lgd1.BoxFace.ColorType = 'truecoloralpha';
lgd1.Box = 'off';
set(gca,'color','white')
% exportgraphics(gca,'R vs T_CEQW241009_CEQW240715_B.tiff','Resolution',600)

%% Rxx vs H

targetTemperatures = [2, 4, 6, 8, 10, 20, 80];
tolerance = 0.1;  % Define a suitable tolerance
T_InSb_Rxx = Data3{:, 1}; 
T_InSbBi_Rxx = Data4{:, 1}; 
T_InSb_Rxy = Data5{:, 1}; 
T_InSbBi_Rxy = Data6{:, 1}; 


% Initialize vectors to hold the corresponding values from columns 2 and 3
B_InSb_Rxx = [];
Rxx_InSb = [];
B_InSbBi_Rxx = [];
Rxx_InSbBi = [];
B_InSb_Rxy = [];
Rxy_InSb = [];
B_InSbBi_Rxy = [];
Rxy_InSbBi = [];

% Loop through the target temperatures and extract corresponding values
for i = 1:length(targetTemperatures)
    % Find the index of the current target temperature
    index_InSb_Rxx = find(abs(T_InSb_Rxx - targetTemperatures(i)) <= tolerance);
    index_InSbBi_Rxx = find(abs(T_InSbBi_Rxx - targetTemperatures(i)) <= tolerance);
    index_InSb_Rxy = find(abs(T_InSb_Rxy - targetTemperatures(i)) <= tolerance);
    index_InSbBi_Rxy = find(abs(T_InSbBi_Rxy - targetTemperatures(i)) <= tolerance);

    if ~isempty(index_InSb_Rxx)
        B_InSb_Rxx{i} = Data3{index_InSb_Rxx, 2};  % Column 2
        Rxx_InSb{i} = Data3{index_InSb_Rxx, 3};  % Column 3
    end

    if ~isempty(index_InSbBi_Rxx)
        B_InSbBi_Rxx{i} = Data4{index_InSbBi_Rxx, 2};  % Column 2
        Rxx_InSbBi{i} = Data4{index_InSbBi_Rxx, 4};  % Column 3
    end

    if ~isempty(index_InSb_Rxy)
        B_InSb_Rxy{i} = Data5{index_InSb_Rxy, 2};  % Column 2
        Rxy_InSb{i} = Data5{index_InSb_Rxy, 3};  % Column 3
    end

    if ~isempty(index_InSbBi_Rxy)
        B_InSbBi_Rxy{i} = Data6{index_InSbBi_Rxy, 2};  % Column 2
        Rxy_InSbBi{i} = Data6{index_InSbBi_Rxy, 4};  % Column 3
    end


end

cmap = jet(length(targetTemperatures)); % 'jet' is a color map that transitions from blue to red
legend_entries = {};

figure4 = figure('WindowState', 'maximized');
for i = 1:length(targetTemperatures)
    plot(B_InSb_Rxx{:,i},Rxx_InSb{:,i},'.-','Color', cmap(i,:),'MarkerSize',14,'LineWidth',2);
    hold on;
    legend_entries{i} = [ num2str(targetTemperatures(i)), ' K']; 
end
hold off;
axis square
xlabel('B (T)')
ylabel('R_{xx} (Ω)') ;
set(gca,'FontSize',22,'FontName','sans serif')
ax = gca;
xlim([-10 10])
ylim([0 max(Rxx_InSb{:,1})+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
lgd1=legend(legend_entries); % Add the legend to the plot in sorted temperature order
set(lgd1,'Position',[0.734071178009941 0.082523377581926 0.139973960878948 0.87931967497387],'EdgeColor',[1 1 1]);
% lgd1.Color = 'none';
% lgd1.BoxFace.ColorType = 'truecoloralpha';
% lgd1.Box = 'off';
set(gca,'color','white')
% exportgraphics(gca,'Rxx vs T_9T_CEQW240715_B.tiff','Resolution',600)


legend_entries = {};
figure5 = figure('WindowState', 'maximized');
for i = 1:length(targetTemperatures)
    plot(B_InSbBi_Rxx{:,i},Rxx_InSbBi{:,i},'.-','Color', cmap(i,:),'MarkerSize',14,'LineWidth',2);
    hold on;
    legend_entries{i} = [ num2str(targetTemperatures(i)), ' K']; 
end
hold off;
axis square
xlabel('B (T)')
ylabel('R_{xx} (Ω)') ;
set(gca,'FontSize',22,'FontName','sans serif')
ax = gca;
xlim([-10 10])
ylim([500 max(Rxx_InSbBi{:,1})+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
lgd1=legend(legend_entries); % Add the legend to the plot in sorted temperature order
set(lgd1,'Position',[0.734071178009941 0.082523377581926 0.139973960878948 0.87931967497387],'EdgeColor',[1 1 1]);
% lgd1.Color = 'none';
% lgd1.BoxFace.ColorType = 'truecoloralpha';
% lgd1.Box = 'off';
set(gca,'color','white')
% exportgraphics(gca,'Rxx vs T_9T_CEQW241009.tiff','Resolution',600)

figure6 = figure('WindowState', 'maximized');
for i = 1:length(targetTemperatures)
    plot(B_InSb_Rxy{:,i},Rxy_InSb{:,i},'.-','Color', cmap(i,:),'MarkerSize',14,'LineWidth',2);
    hold on;
    legend_entries{i} = [ num2str(targetTemperatures(i)), ' K']; 
end
hold off;
axis square
xlabel('B (T)')
ylabel('R_{xy} (Ω)') ;
set(gca,'FontSize',22,'FontName','sans serif')
ax = gca;
xlim([-10 10])
ylim([min(Rxy_InSb{:,1})-100 max(Rxy_InSb{:,1})+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
lgd1=legend(legend_entries); % Add the legend to the plot in sorted temperature order
set(lgd1,'Position',[0.734071178009941 0.082523377581926 0.139973960878948 0.87931967497387],'EdgeColor',[1 1 1]);
% lgd1.Color = 'none';
% lgd1.BoxFace.ColorType = 'truecoloralpha';
% lgd1.Box = 'off';
set(gca,'color','white')
% exportgraphics(gca,'Rxy vs T_9T_CEQW240715_B.tiff','Resolution',600)


legend_entries = {};
figure7 = figure('WindowState', 'maximized');
for i = 1:length(targetTemperatures)
    plot(B_InSbBi_Rxy{:,i},-1*Rxy_InSbBi{:,i},'.-','Color', cmap(i,:),'MarkerSize',14,'LineWidth',2);
    hold on;
    legend_entries{i} = [ num2str(targetTemperatures(i)), ' K']; 
end
hold off;
axis square
xlabel('B (T)')
ylabel('R_{xy} (Ω)') ;
set(gca,'FontSize',22,'FontName','sans serif')
ax = gca;
xlim([-10 10])
ylim([min(Rxy_InSbBi{:,1})-100 max(Rxy_InSbBi{:,1})+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
lgd1=legend(legend_entries); % Add the legend to the plot in sorted temperature order
set(lgd1,'Position',[0.734071178009941 0.082523377581926 0.139973960878948 0.87931967497387],'EdgeColor',[1 1 1]);
% lgd1.Color = 'none';
% lgd1.BoxFace.ColorType = 'truecoloralpha';
% lgd1.Box = 'off';
set(gca,'color','white')
% exportgraphics(gca,'Rxy vs T_9T_CEQW241009.tiff','Resolution',600)

