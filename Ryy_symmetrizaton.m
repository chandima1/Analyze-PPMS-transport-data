clc
clear all
close all

% Setup import options for selective columns
opts = delimitedTextImportOptions("NumVariables", 42);

% Specify range, delimiter, and file-level properties
opts.DataLines = [32, Inf];  % Adjust the start row if needed
opts.Delimiter = ",";  % Change this if your file uses a different delimiter like tabs '\t'


%% Ryy symmetrization

% Define the file path
filePath1 = 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\Chan2_CEQW241010_Ryy_vs_H_9T__5uA_1.8Kto80K_Chan3_CEQW241009_Ryy_vs_H_9T__5uA_1.8Kto80K_new_00001.dat';
% filePath2 = 'E:\Transport\CEIQW240911_B and CEQW240715_B\Chan3_Ryy_vs_H_complete_9T__5uA_2Kto50K_InSbBi240911_B_Chan2_not working_InSb240715_B_00001.dat';

% Import the data
data1 = readtable(filePath1, opts);
% data2 = readtable(filePath2, opts);

% Clear temporary variables
clear opts;

T1 = str2double(data1{:, 4});
B1 = str2double(data1{:, 5})/10000;
Res2 = str2double(data1{:,21});
Res3 = str2double(data1{:,22});


%%

% Initialize variables
threshold = 0.15; % 0.5 °C threshold
temperature_blocks = {}; % Cell array to hold the blocks
current_block = []; % Array to hold the current block of temperatures


% Iterate through the temperature data
for i = 1:length(T1)
    if isempty(current_block)
        % Start a new block if current block is empty
        current_block = T1(i);
    else
        % Check the difference with the last temperature in the current block
        if abs(T1(i) - current_block(end)) <= threshold
            % If within threshold, add to current block
            current_block = [current_block, T1(i)];
        else
            % If out of threshold, store the current block and start a new one
            temperature_blocks{end + 1} = current_block; % Save the block
            current_block = T1(i); % Start a new block
        end
    end
end

% Save the last block if not already saved
if ~isempty(current_block)
    temperature_blocks{end + 1} = current_block; % Save the final block
end

No_of_temp=length(temperature_blocks);
no_of_blocks = cellfun(@length, temperature_blocks);
iii=0;
Sym_Rxx2 = cell(1, length(no_of_blocks)); % Cell array for symmetrized Rxx
Sym_Rxx3 = cell(1, length(no_of_blocks)); % Cell array for symmetrized Rxx
Sym_B = cell(1, length(no_of_blocks));   % Cell array for symmetrized B

for i=1:length(no_of_blocks)
    current_block_length = length(temperature_blocks{i}); % Length of the current block
    half_length = floor(current_block_length / 2); % Calculate half the length for symmetry

    % Preallocate symmetrized arrays for the current block
    Sym_Rxx2{i} = zeros(1, half_length);
    Sym_B{i} = zeros(1, half_length);
    for ii=1:half_length
        Sym_Rxx2{i}(ii)=(Res2(ii+iii)+Res2(length(temperature_blocks{i})/2+ii+iii-1))/2;
        Sym_Rxx3{i}(ii)=(Res3(ii+iii)+Res3(length(temperature_blocks{i})/2+ii+iii-1))/2;
        % Sym_B{i}(ii)=(B(ii)+B(length(temperature_blocks{i})-ii+1))/2;
        % Sym_B{i}(ii)=B(ii+iii);
        if B1(ii+iii)>0
            Sym_B{i}(ii)=(B1(ii+iii)+abs(B1(length(temperature_blocks{i})/2+ii+iii-1)))/2;
        else
            Sym_B{i}(ii)=-1*(abs(B1(ii+iii))+B1(length(temperature_blocks{i})/2+ii+iii-1))/2;
        end
    end
    iii=iii+length(temperature_blocks{i});
end

%%
figure1=figure('WindowState','maximized');
hold on; % Hold on to add multiple plots to the same figure
legend_entries = cell(1, No_of_temp); 
temperatures = zeros(1, No_of_temp);
for j = 1:No_of_temp
    temperatures(j) = round(temperature_blocks{j}(1), 2); % Assuming temperature is the same for the entire block
end

% Sort temperatures in ascending order and get the sorting index
[sorted_temperatures, sort_idx] = sort(temperatures);

% Get a colormap for the plots, where each row is an RGB triplet
cmap = jet(No_of_temp); % 'jet' is a color map that transitions from blue to red

% Loop through each block for plotting based on sorted temperature index
for j = 1:No_of_temp
    % Get the sorted index
    idx = sort_idx(j);
    plot(Sym_B{idx}(1:180), Sym_Rxx2{idx}(1:180), '.', 'MarkerSize', 14, 'Color', cmap(j,:));
    legend_entries{j} = ['T = ', num2str(sorted_temperatures(j)), ' K']; % Assign legend entry with sorted temperature
end

axis square; % Set the axes to be square
xlabel('Magnetic Field (T)'); % Add x-axis label
ylabel('Symmetrized R_{yy} (Ohms)');   % Add y-axis label
lgd1=legend(legend_entries); % Add the legend to the plot in sorted temperature order
set(lgd1,'Position',[0.734071178009941 0.082523377581926 0.139973960878948 0.87931967497387]);
% lgd.Color = 'none';
% lgd.BoxFace.ColorType = 'truecoloralpha';
% lgd.Box = 'off';
hold off; % Release the hold
set(gca,'FontSize',22,'FontName','sans serif');
axis square
ax = gca;
xlim([-10 10])
% ylim([min(Sym_Rxx2{No_of_temp}(:))-100 max(Sym_Rxx2{1}(:))+100])
ylim([-2.1e5 1.4e6])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
set(gca,'color','white')
exportgraphics(gca,'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\CEQW241010_Symmetrized_Ryy.tiff','Resolution',600)





figure2=figure('WindowState','maximized');
hold on; % Hold on to add multiple plots to the same figure
legend_entries = cell(1, No_of_temp); 
temperatures = zeros(1, No_of_temp);
for j = 1:No_of_temp
    temperatures(j) = round(temperature_blocks{j}(1), 2); % Assuming temperature is the same for the entire block
end
[sorted_temperatures, sort_idx] = sort(temperatures);
cmap = jet(No_of_temp); 
for j = 1:No_of_temp 
    % Get the sorted index
    if j == 9
        continue;
    end

    idx = sort_idx(j);
    plot(Sym_B{idx}(1:180), Sym_Rxx3{idx}(1:180), '.', 'MarkerSize', 14, 'Color', cmap(j,:));
    legend_entries{j} = ['T = ', num2str(sorted_temperatures(j)), ' K']; % Assign legend entry with sorted temperature
    % pause(2)
end
axis square; % Set the axes to be square
xlabel('Magnetic Field (T)'); % Add x-axis label
ylabel('Symmetrized R_{yy} (Ohms)');   % Add y-axis label

% Remove empty entries from legend
legend_entries = legend_entries(~cellfun('isempty', legend_entries));
lgd2 = legend(legend_entries); % Create legend with non-empty entries
set(lgd2,'Position',[0.734071178009941 0.082523377581926 0.139973960878948 0.87931967497387]);
% lgd.Color = 'none';
% lgd.BoxFace.ColorType = 'truecoloralpha';
% lgd.Box = 'off';
hold off; % Release the hold
set(gca,'FontSize',22,'FontName','sans serif')
axis square
ax = gca;
xlim([-10 10])
ylim([min(Sym_Rxx3{No_of_temp}(:))-100 max(Sym_Rxx3{1}(:))+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
set(gca,'color','white')
exportgraphics(gca,'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\CEQW241009_Symmetrized_Ryy.tiff','Resolution',600)



%%

% Initialize a cell array to hold the combined data
combined_data = cell(sum(cellfun(@length, temperature_blocks)), 4); % 3 columns: Temperature, Magnetic Field, Resistance

row = 1; % Row index for combined_data

for j = 1:No_of_temp
    idx = sort_idx(j); % Get the sorted index
    T = round(temperature_blocks{idx}(1), 3); % Temperature for the current block
    
    % Ensure we're using the full length of Sym_B and Sym_Rxx2 for the current block
    for k = 1:length(Sym_B{idx}) % Loop through the magnetic field data
        % Check if the index k is within the range
        if k <= length(Sym_B{idx}) && k <= length(Sym_Rxx3{idx}) && k <= length(Sym_Rxx2{idx})
            combined_data(row, :) = {T, Sym_B{idx}(k), Sym_Rxx2{idx}(k), Sym_Rxx3{idx}(k)}; % Fill in the temperature, B, and resistance
            row = row + 1; % Increment the row index
        end
    end
end

% Convert to a table for better CSV formatting
combined_table = cell2table(combined_data, 'VariableNames', {'Temperature', 'MagneticField', 'Symmetrized_Ryy_Chan2', 'Symmetrized_Ryy_Chan3'});

% Write to CSV
writetable(combined_table, 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\Symmetrized_Ryy_9T_Chan2_CEQW241010_Chan3_CEQW241009.csv');

