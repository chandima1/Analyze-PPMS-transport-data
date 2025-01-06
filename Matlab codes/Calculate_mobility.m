clc;
clear all;
close all;

%% Set up the Import Options and import the data

filename1 = 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\Symmetrized_Rxx_9T_Chan2_CEQW241010_Chan3_CEQW241009.csv';
filename2 = 'E:\Transport\Chan2_CEQW241010_Chan3_CEQW241009\Symmetrized_Ryy_9T_Chan2_CEQW241010_Chan3_CEQW241009.csv';

Data1 = readtable(filename1);
Data2 = readtable(filename2);


%% Clear temporary variables
T1=Data1{:,1};
% T2=Data2{:,1};
B1=Data1{:,2};
B2=Data2{:,2};
Rxx_Chan2=Data1{:,3};
Ryy_Chan2=Data2{:,3};
Rxx_Chan3=Data1{:,4};
Ryy_Chan3=Data2{:,4};
upper_bound=2.1;
indices = find(T1 < upper_bound);
% indices_Chan3 = find(T2 < upper_bound);

Sym_B_Chan2=Data1{min(indices):max(indices),2};
Sym_B_Chan3=Data1{min(indices):max(indices),2};
Sym_Rxx_Chan2=Data1{min(indices):max(indices),3};
Sym_Rxx_Chan3=Data1{min(indices):max(indices),4};
Sym_Ryy_Chan2=Data2{min(indices):max(indices),3};
Sym_Ryy_Chan3=Data2{min(indices):max(indices),4};

% Example sets of known values
R_vertical_values_Chan2 = Sym_Ryy_Chan2;
R_horizontal_values_Chan2 = Sym_Rxx_Chan2;

R_vertical_values_Chan3 = Sym_Ryy_Chan3;
R_horizontal_values_Chan3 = Sym_Rxx_Chan3;

% Preallocate array for R_S results
R_S_results_Chan2 = zeros(size(R_vertical_values_Chan2));

R_S_results_Chan3 = zeros(size(R_vertical_values_Chan3));

% Loop through each pair of values
for i = 1:length(R_vertical_values_Chan2)
    R_vertical_Chan2 = R_vertical_values_Chan2(i);
    R_horizontal_Chan2 = R_horizontal_values_Chan2(i);
    
    % Initial guess for Rs
    initial_guess = 10000;

    % Use fsolve to find the root
    options = optimset('Display', 'off'); % Suppress output
    R_S_results_Chan2(i) = fsolve(@(RS) equations(RS, R_vertical_Chan2, R_horizontal_Chan2), initial_guess, options);
end


plot(Sym_B_Chan2,R_S_results_Chan2,'r.','MarkerSize',15)
hold on
plot(Sym_B_Chan2,Sym_Rxx_Chan2,'g.','MarkerSize',15)
plot(Sym_B_Chan2,Sym_Ryy_Chan2,'b.','MarkerSize',15)
hold off
legend('Sheet resistance','Symmetrized Rxx','Symmetrized Ryy')
[~, closest_index]=(min(abs((Sym_B_Chan2))));
set(gca,'FontSize',22,'FontName','Sans Serif')

q=1.602e-19;
pre_factor=1e4; % to change to cm^2
n= -4.4396e+15 % find from the Rxy_symmetrization.m code. Here, n is in m-2
Sheet_res=R_S_results_Chan2(closest_index)
Mobility_in_m2_Chan2=1/(q*n*Sheet_res) % In m^2/V.s
Mobility_in_cm2_CHan2=1e4/(q*n*Sheet_res) % In cm^2/V.s


for i = 1:length(R_vertical_values_Chan3)
    R_vertical_Chan3 = R_vertical_values_Chan3(i);
    R_horizontal_Chan3 = R_horizontal_values_Chan3(i);
    
    % Initial guess for Rs
    initial_guess = 10000;

    % Use fsolve to find the root
    options = optimset('Display', 'off'); % Suppress output
    R_S_results_Chan3(i) = fsolve(@(RS) equations(RS, R_vertical_Chan3, R_horizontal_Chan3), initial_guess, options);
end

figure
plot(Sym_B_Chan3,R_S_results_Chan3,'r.','MarkerSize',15)
hold on
plot(Sym_B_Chan3,Sym_Rxx_Chan3,'g.','MarkerSize',15)
plot(Sym_B_Chan3,Sym_Ryy_Chan3,'b.','MarkerSize',15)
hold off
legend('Sheet resistance','Symmetrized Rxx','Symmetrized Ryy')
[~, closest_index]=(min(abs((Sym_B_Chan2))));
set(gca,'FontSize',22,'FontName','Sans Serif')

q=1.602e-19;
pre_factor=1e4; % to change to cm^2
n_Chan3= -4.4396e+15 % find from the Rxy_symmetrization.m code. Here, n is in m-2
Sheet_res=R_S_results_Chan3(closest_index)
Mobility_in_m2_Chan3=1/(q*n_Chan3*Sheet_res) % In m^2/V.s
Mobility_in_cm2_Chan3=1e4/(q*n_Chan3*Sheet_res) % In cm^2/V.s


% Define the function to solve
function F = equations(RS, R_vertical, R_horizontal)
    F = exp(-pi * R_vertical / RS) + exp(-pi * R_horizontal / RS) - 1;
end
