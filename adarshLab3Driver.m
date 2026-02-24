clc;
clear;
close all;

%% RW Calc

data = readmatrix("2026_02_10_002_04_RW_T20t210");

time = data(:, 1) / 1000; % [s]
comTorque = data(:, 2) / 1000;
angVel = data(:, 3) .* (2 * pi / 60);
current = data(:, 4);


Kt = 33.5 / 1000; % from data sheet

range = (time >= 1.2) & (time <= 3.8); % by inspection
t_lin = time(range);
angVel_lin = angVel(range);
current_lin = current(range);

torque_app = mean(current_lin) * Kt;
p = polyfit(t_lin, angVel_lin, 1);
angAcc = p(1);

I_rw = torque_app / angAcc;

%% SC Gains Calc

data = readmatrix("2026_02_24_002_04_SC_T10t5");


time_sc = data(:, 1) / 1000; % [s]
comTorque_sc = data(:, 2) / 1000;
angVel_sc = data(:, 3) .* (2 * pi / 60);
current_sc = data(:, 4);

Kt_sc = 25.5 / 1000; % from data sheet

range_sc = (time_sc >= 2) & (time_sc <= 5); % by inspection
t_lin_sc = time_sc(range_sc);
angVel_lin_sc = angVel_sc(range_sc);
current_lin_sc = current_sc(range_sc);

torque_app_sc = mean(current_lin_sc) * Kt_sc;
p = polyfit(t_lin_sc, angVel_lin_sc, 1);
angAcc_sc = p(1);

I_sc = abs(torque_app_sc / angAcc_sc);

t_max = 1.5; % [s]
Mp_max = 0.10;

zeta_min = sqrt((log(Mp_max)^2) / (pi^2 + log(Mp_max)^2));
zeta = 0.707; 

wn_min = 3 / (zeta * t_max);
wn = 3; % select wn > wn_min

K1 = I_sc * (wn^2); % [Nm/rad]
K2 = I_sc * (2 * zeta * wn); % [Nm/(rad/s)]

K1 = K1 * 1000; % [mNm/rad]
K2 = K2 * 1000; % [mNm/(rad/s)]

%% Gains Analysis on Chosen / Tests

folderPath = 'C:\Users\adbo5750exc\iCloudDrive\ASEN 3801\Lab3\Gains Testing'; % NOTE. THIS WILL BE CHANGED PER PERSON WHEN RUNNING THIS CODE
files = dir(folderPath);

dataStr = struct();

l = 0;

for k = 1:length(files)

    if files(k).isdir
        continue;
    end

    l = l + 1; 

    baseFileName = files(k).name;
    fullFileName = fullfile(folderPath, baseFileName);

    data_temp = readmatrix(fullFileName);
    first_real_idx = find(data_temp(:, 1) > 0, 1);

    if ~isempty(first_real_idx)
        start_time = data_temp(first_real_idx, 1);
        data_temp(first_real_idx:end, 1) = data_temp(first_real_idx:end, 1) - start_time;
    end

    data_temp = extract_full_period(data_temp);

    dataStr(l).time = (data_temp(:, 1) / 1000) ;
    dataStr(l).ref_pos = max(data_temp(:, 2));
    dataStr(l).measured_pos = data_temp(:, 3);
    dataStr(l).current = data_temp(:, 4);
    dataStr(l).Kp = data_temp(end, 5);
    dataStr(l).Kd = data_temp(end, 6);
end


T = struct2table(dataStr);
T_sorted = sortrows(T, {'Kp', 'Kd'});
finDataStr = table2struct(T_sorted);
finDataStr = finDataStr';

for i = 1:length(finDataStr)
    time = finDataStr(i).time;
    measuredPos = finDataStr(i).measured_pos;
    refPos = finDataStr(i).ref_pos;
    dispText = "Kp = " + num2str(finDataStr(i).Kp) + " Kd = " + num2str(finDataStr(i).Kd);
    
    figure(100 + i)
    plot(time, measuredPos, 'b', time, refPos * ones(size(time)), 'r--');
    text(0, min(measuredPos), dispText)
    xlabel('Time (s)');
    ylabel('Position (rad)');
    title(['Position Response for Test ' num2str(i)]);
    legend('Measured Position', 'Reference Position');
    grid on;
end

save_figs = 0; % Set to 1 to save, 0 to skip

if save_figs == 1
    figs = findobj('Type', 'figure'); 
    
    for i = 1:length(figs)
        % Saves each figure as a 300 DPI PNG using its figure number
        filename = sprintf('Figure_%d.png', figs(i).Number);
        exportgraphics(figs(i), filename, 'Resolution', 300);
    end
end

%% function to extract 1 period. DID NOT WRITE

function period_data = extract_full_period(data)
    % EXTRACT_FULL_PERIOD Cleans the matrix to include exactly one full 
    % period of oscillation. It uses zero-crossings of the Reference 
    % Position (Column 2) to determine the cycle boundaries.
    % 
    % Usage: 
    %   cleaned_matrix = extract_full_period(my_data_matrix);

    % 1. Remove the initialization dummy row (where time == 0.000)
    first_real_idx = find(data(:, 1) > 0, 1);
    if isempty(first_real_idx)
        error('No valid timestamp found in the first column.');
    end
    clean_data = data(first_real_idx:end, :);
    
    % 2. Extract the reference signal (Column 2)
    ref_signal = clean_data(:, 2);
    
    % 3. Mean-center the signal to find crossings easily
    % By subtracting the mean, the signal will alternate between positive 
    % and negative values, making transitions explicitly clear.
    ref_centered = ref_signal - mean(ref_signal);
    
    % 4. Find indices where the sign changes (the zero crossings)
    % A non-zero difference means a transition happened between index i and i+1
    crossings = find(diff(sign(ref_centered)) ~= 0);
    
    % 5. Isolate one full cycle (two consecutive half-cycles)
    % We discard the data before the 1st crossing, as it is usually just a 
    % partial cycle or a starting steady-state hold.
    if length(crossings) >= 3
        % Start at the beginning of the first full half-cycle
        start_idx = crossings(1) + 1;
        
        % End exactly at the end of the second full half-cycle
        end_idx = crossings(3);
        
        period_data = clean_data(start_idx:end_idx, :);
        
    elseif length(crossings) == 2
        warning('Only one full half-cycle found. Returning the rest of the data.');
        start_idx = crossings(1) + 1;
        period_data = clean_data(start_idx:end, :);
        
    else
        warning('No complete oscillation cycles detected. Returning all normalized data.');
        period_data = clean_data;
    end
    
    % 6. Normalize the time column so the extracted period starts exactly at t = 0
    period_data(:, 1) = period_data(:, 1) - period_data(1, 1);
end
