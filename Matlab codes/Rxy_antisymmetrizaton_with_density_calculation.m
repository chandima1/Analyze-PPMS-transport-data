clc
clear all
close all

% Setup import options for selective columns
opts = delimitedTextImportOptions("NumVariables", 42);

% Specify range, delimiter, and file-level properties
opts.DataLines = [32, Inf];  % Adjust the start row if needed
opts.Delimiter = ",";  % Change this if your file uses a different delimiter like tabs '\t'


%% Rxx symmetrization

% Define the file path
filePath1 = 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\Chan2_CEQW241010_Rxy_vs_H_9T__5uA_1.8Kto80K_Chan3_CEQW241009_Rxy_vs_H_9T__5uA_1.8Kto80K_new_00001.dat';


% Import the data
data1 = readtable(filePath1, opts);
% data2 = readtable(filePath2, opts);

% Clear temporary variables
clear opts;

T1 = str2double(data1{:, 4});
B1 = str2double(data1{:, 5})/10000;
Res2 = str2double(data1{:,21});
Res3 = -1*str2double(data1{:,22});


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
Sym_Rxy2 = cell(1, length(no_of_blocks)); % Cell array for symmetrized Rxx
Sym_Rxy3 = cell(1, length(no_of_blocks)); % Cell array for symmetrized Rxx
Sym_B = cell(1, length(no_of_blocks));   % Cell array for symmetrized B

for i=1:length(no_of_blocks)
    current_block_length = length(temperature_blocks{i}); % Length of the current block
    half_length = floor(current_block_length / 2); % Calculate half the length for symmetry

    % Preallocate symmetrized arrays for the current block
    Sym_Rxy2{i} = zeros(1, half_length);
    Sym_B{i} = zeros(1, half_length);
    for ii=1:half_length
        Sym_Rxy2{i}(ii)=(Res2(ii+iii)-Res2(current_block_length/2+ii+iii-1))/2;
        Sym_Rxy3{i}(ii)=(Res3(ii+iii)-Res3(current_block_length/2+ii+iii-1))/2;
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
    plot(Sym_B{idx}(1:180), Sym_Rxy2{idx}(1:180), '.', 'MarkerSize', 14, 'Color', cmap(j,:));
    legend_entries{j} = ['T = ', num2str(sorted_temperatures(j)), ' K']; % Assign legend entry with sorted temperature
end

axis square; % Set the axes to be square
xlabel('Magnetic Field (T)'); % Add x-axis label
ylabel('Antisymmetrized R_{xy} (Ohms)');   % Add y-axis label
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
% ylim([min(Sym_Rxy2{1}(:))-100 max(Sym_Rxy2{1}(:))+100])
ylim([-11000 11000])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
set(gca,'color','white')
% exportgraphics(gca,'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\CEQW241010_Antisymmetrized_Rxy.tiff','Resolution',600)

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
    idx = sort_idx(j);
    plot(Sym_B{idx}(1:180), Sym_Rxy3{idx}(1:180), '.', 'MarkerSize', 14, 'Color', cmap(j,:));
    legend_entries{j} = ['T = ', num2str(sorted_temperatures(j)), ' K']; % Assign legend entry with sorted temperature
end
axis square; % Set the axes to be square
xlabel('Magnetic Field (T)'); % Add x-axis label
ylabel('Antisymmetrized R_{xy} (Ohms)');   % Add y-axis label
lgd2=legend(legend_entries); % Add the legend to the plot in sorted temperature order
set(lgd2,'Position',[0.734071178009941 0.082523377581926 0.139973960878948 0.87931967497387]);
% lgd.Color = 'none';
% lgd.BoxFace.ColorType = 'truecoloralpha';
% lgd.Box = 'off';
hold off; % Release the hold
set(gca,'FontSize',22,'FontName','sans serif')
axis square
ax = gca;
xlim([-10 10])
ylim([min(Sym_Rxy3{1}(:))-100 max(Sym_Rxy3{1}(:))+100])
currentXLim = xlim;
currentYLim = ylim;
ax=gca;
ax.Box = 'off';
ax.LineWidth = 1.5; % Make border lines thicker
xline(currentXLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
yline(currentYLim(2),'k','LineWidth',2,'HandleVisibility', 'off')
set(gca,'color','white')
% exportgraphics(gca,'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\CEQW241009_Antisymmetrized_Rxy.tiff','Resolution',600)

%% Find the density from the curve
Rxy3=Sym_Rxy3{1}(:);
Rxy2=Sym_Rxy2{1}(:);
B=Sym_B{1}(:);

% Fit a linear model
p_chan3 = polyfit(B, Rxy3, 1); % p(1) = slope, p(2) = intercept
p_chan2 = polyfit(B, Rxy2, 1); % p(1) = slope, p(2) = intercept

% Generate fitted values
Rxy3_fit = polyval(p_chan3, B);
Rxy2_fit = polyval(p_chan2, B);


figure 
plot(B,Rxy2,'b.')
hold on
plot(B, Rxy2_fit, 'r-', 'LineWidth', 1, 'DisplayName', 'Fitted Line');
xlabel('Magnetic Field (B)');
ylabel('R_{xy}');
title('Linear Fit to Data');
legend show;
grid on;
hold off;
axis square;
t=20e-9%thickness in meter;
grad=p_chan3(1);
Carrier_density_3D=1/(grad*1.602e-19*t); % In m-3 
% the 3D to 2D conversion => n(in m^-2)=n(in m^-3) * t(in m)
Carrier_density_2D_in_m2_Chan2=Carrier_density_3D*t % In m-2
Carrier_density_2D_in_cm2_Chan2=Carrier_density_3D*t*1e-4 % In cm-2


figure 
plot(B,Rxy3,'b.')
hold on
plot(B, Rxy3_fit, 'r-', 'LineWidth', 1, 'DisplayName', 'Fitted Line');
xlabel('Magnetic Field (B)');
ylabel('R_{xy}');
title('Linear Fit to Data');
legend show;
grid on;
hold off;
axis square;
t=20e-9%thickness in meter;
grad=p_chan3(1);
Carrier_density_3D=1/(grad*1.602e-19*t); % In m-3 
% the 3D to 2D conversion => n(in m^-2)=n(in m^-3) * t(in m)
Carrier_density_2D_in_m2_Chan3=Carrier_density_3D*t % In m-2
Carrier_density_2D_in_cm2_Chan3=Carrier_density_3D*t*1e-4 % In cm-2



%% Fit 2 band model
% Two-band model function
% Rxy = (-B/e)*((n1*u1^2+n2*u2^2+B^2*u1^2*u2^2*(n1+n2))/((n1*u1+n2u2)^2+B^2*u1^2*u2^2*(n1+n2)^2));

Rxy2=Sym_Rxy2{1}(:);
B=Sym_B{1}(:);

two_band_model = @(params, B) ...
    (-B / params(1)) .* ...
    ((params(2) * params(3)^2 + params(4) * params(5)^2 + ...
      B.^2 * params(3)^2 * params(5)^2 * (params(2) + params(4))) ./ ...
     ((params(2) * params(3) + params(4) * params(5))^2 + ...
      B.^2 * params(3)^2 * params(5)^2 * (params(2) + params(4))^2));

% Initial parameter guess: [e, n1, u1, n2, u2]
initial_params = [1.6e-19, 2.1e15, 0.597 ,-1.07e15, -0.279]; % Example values

opts = optimset('Display', 'off');
[param_fit, resnorm] = lsqcurvefit(two_band_model, initial_params, B, Rxy2, [], [], opts);


% Calculate fitted Rxy using the fitted parameters
Rxy2_fit = two_band_model(param_fit, B);

% Plot the results
figure;
plot(B, Rxy2, 'ro', 'MarkerFaceColor', 'r'); % Original data
hold on;
plot(B, Rxy2_fit, 'b-'); % Fitted curve
xlabel('Magnetic Field (B) [T]');
ylabel('Hall Resistance (Rxy) [Ω]');
title('Two-Band Model Fit to Rxy vs. B Data');
legend('Data', 'Fit');
grid on;
hold off;

% fit linear for a small range
upper_bound=1; % in Tesla
indices = find(B > -1*upper_bound & B < upper_bound);
B_short_range=B(min(indices):max(indices));
Rxy_short_range=Rxy2(min(indices):max(indices));
% Fit a linear model
p_short_range = polyfit(B_short_range, Rxy_short_range, 1); % p(1) = slope, p(2) = intercept

% Generate fitted values
Rxy2_fit_short_range = polyval(p_short_range, B_short_range);

figure
plot(B_short_range,Rxy_short_range,'b.')
hold on
plot(B_short_range,Rxy2_fit_short_range,'r-.')
hold off
axis square
grad_short_range=p_short_range(1);
Carrier_density_3D_short_range=1/(grad_short_range*1.602e-19*t); % In m-3 
% the 3D to 2D conversion => n(in m^-2)=n(in m^-3) * t(in m)
Carrier_density_2D_in_m2_short_range=Carrier_density_3D_short_range*t % In m-2
Carrier_density_2D_in_cm2_short_range=Carrier_density_3D_short_range*t*1e-4 % In cm-2

%% 


% Initialize a cell array to hold the combined data
combined_data = cell(sum(cellfun(@length, temperature_blocks)), 4); % 3 columns: Temperature, Magnetic Field, Resistance

row = 1; % Row index for combined_data

for j = 1:No_of_temp
    idx = sort_idx(j); % Get the sorted index
    T = round(temperature_blocks{idx}(1), 2); % Temperature for the current block
    
    % Ensure we're using the full length of Sym_B and Sym_Rxx2 for the current block
    for k = 1:length(Sym_B{idx}) % Loop through the magnetic field data
        % Check if the index k is within the range
        if k <= length(Sym_B{idx}) && k <= length(Sym_Rxy2{idx})
            combined_data(row, :) = {T, Sym_B{idx}(k), Sym_Rxy2{idx}(k),Sym_Rxy3{idx}(k)}; % Fill in the temperature, B, and resistance
            row = row + 1; % Increment the row index
        end
    end
end

% Convert to a table for better CSV formatting
combined_table = cell2table(combined_data, 'VariableNames', {'Temperature', 'MagneticField', 'Symmetrized_Rxy_Chan2', 'Symmetrized_Rxy_Chan3'});

% Write to CSV
writetable(combined_table, 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\Symmetrized_Rxy_9T_Chan2_CEQW241010_Chan3_CEQW241009.csv');
