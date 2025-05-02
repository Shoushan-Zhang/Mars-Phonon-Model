function [depth_min,depth_max,rms_perturbation,scale_length,P_velocity,S_velocity,Qalpha_depth_min,Qalpha_depth_max,Qalpha]=funPlotDoFile(filename)
    % funPlotDoFile extracts and plots scattering parameters and Qalpha from a given file
    %        Also it plots the uncertainty of RMS and Qi
    %
    % Inputs:
    %   filename - The path to the configuration file (e.g., 'do.photon')
    %
    % Example:
    %   plot_scattering_parameters_with_Qalpha('do.photon')

    % Open the file

Qi=[10000;10000;2500;660;1520;4280;8100];
Qi= 9/4 * Qi * 3.5;


eps = ...
[
    0.0844
    0.0194
    0.0090
    0.0098
    0.0082
    0.0082
    0.0082];

depth_index = ...
[
10
70
150
300
600
1050
1500];

    fid = fopen(filename, 'r');
    
    % Check if the file opened successfully
    if fid == -1
        error('Could not open the file. Please check the filename and path.');
    end
    
    % Initialize arrays to store the parameters
    depth_min = [];
    depth_max = [];
    rms_perturbation = [];
    scale_length = [];
    P_velocity = [];
    S_velocity = [];
    Qalpha_depth_min = [];
    Qalpha_depth_max = [];
    Qalpha = [];
    
    % Initialize line counter
    line_num = 0;
    
    % Read through the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        line_num = line_num + 1;
        
        % Display the current line number and line content for debugging
        % disp(['Reading line ', num2str(line_num), ': ', line]);

        
        
        % Process scattering parameters
        if contains(line, 'min scat depth')

            % Extract numerical value before any comment
            numeric_part = strtok(line, '!');
            numeric_value = str2double(strtrim(numeric_part));

            % disp('Found scattering parameters...');
            depth_min(end+1) = numeric_value;
            depth_max(end+1) = str2double(strtrim(strtok(fgetl(fid), '!')));
            fgetl(fid); % Skip 'max scat range from source'
            velocities = str2num(strtok(fgetl(fid), '!')); %#ok<ST2NM> % Read reference P & S velocities
            P_velocity(end+1) = velocities(1);
            S_velocity(end+1) = velocities(2);
            fgetl(fid); % Skip 'relative size of density perturbation'
            rms_perturbation(end+1) = str2double(strtrim(strtok(fgetl(fid), '!')));
            scale_length(end+1) = str2double(strtrim(strtok(fgetl(fid), '!')));
        
        % Process Qalpha parameters
        elseif contains(line, 'depth of Q layer')
            % disp('Found Qalpha layer...');

            line = strtok(line,'!');
            numeric_parts = strsplit(strtrim(line));
            numbers = str2double(numeric_parts);
            Qalpha_depth_min(end+1) = numbers(1);
            Qalpha_depth_max(end+1) = numbers(2);
            

        elseif contains(line, 'Qalpha')
            % disp('Found Qalpha value')

            numeric_part = strtok(line, '!');
            numeric_value = str2double(strtrim(numeric_part));
            Qalpha(end+1) = numeric_value;

        end
    end
    
    % Close the file
    fclose(fid);
    
    % Check if any data was read
    if isempty(depth_min)
        error('No scattering parameters were found in the file.');
    end
    if isempty(Qalpha)
        warning('No Qalpha parameters were found in the file.');
    end
    
    % Create a figure with subplots
   figure;
    
    % Subplot 1: RMS Perturbation vs. Depth
    subplot(1, 4, 1); % 2x2 grid, 1st subplot
   
    for i = 1:length(depth_min)
        plot([rms_perturbation(i), rms_perturbation(i)], [depth_min(i), depth_max(i)],  'LineWidth', 3,'Marker','o'); hold on;
    end
    plot(0.8 * eps,depth_index,'LineStyle','--','Color',[0.5, 0.5, 0.5],'LineWidth',3)
    plot(1.2 * eps,depth_index,'LineStyle','--','Color',[0.5, 0.5, 0.5],'LineWidth',3)
    ylabel('Depth (km)');
    xlabel('RMS');
    title('RMS');
    grid on;
    set(gca,'YDir','reverse','FontSize',16,'FontWeight','bold')
    ylim([-10 1600])
    xlim("auto")
    set(gca,"XMinorTick","on")
    set(gca,"YMinorTick","on")
    set(gca,"XMinorGrid","on")
    set(gca,"YMinorGrid","on")
    

    % Subplot 2: Scale Length vs. Depth
    subplot(1, 4, 2); % 2x2 grid, 2nd subplot

    for i = 1:length(depth_min)
        plot([scale_length(i), scale_length(i)], [depth_min(i), depth_max(i)],  'LineWidth', 3,'Marker','o'); hold on;
    end
    ylabel('Depth (km)');
    xlabel('Scale Length (km)');
    title('Scale Length');
    grid on;
    set(gca,'YDir','reverse','FontSize',16,'FontWeight','bold')
    ylim([-10 1600])
    xlim([min(scale_length)-0.1 max(scale_length)+0.1])
    set(gca,"XMinorTick","on")
    set(gca,"YMinorTick","on")
    set(gca,"XMinorGrid","on")
    set(gca,"YMinorGrid","on")

    % Subplot 3: Reference P & S Velocities vs. Depth
    subplot(1, 4, 3); % 2x2 grid, 3rd subplot

    for i = 1:length(depth_min)
        plot([P_velocity(i), P_velocity(i)], [depth_min(i), depth_max(i)],  'r', 'LineWidth', 3,'Marker','o'); hold on;
        plot([S_velocity(i), S_velocity(i)], [depth_min(i), depth_max(i)],  'b', 'LineWidth', 3,'Marker','o');
    end
    ylabel('Depth (km)');
    xlabel('Velocity (km/s)');
    title('Reference Velocities');
    grid on;
    set(gca,'YDir','reverse','FontSize',16,'FontWeight','bold')
    ylim([-10 1600])
    xlim([min(S_velocity)-1 max(P_velocity)+1])
    set(gca,"XMinorTick","on")
    set(gca,"YMinorTick","on")
    set(gca,"XMinorGrid","on")
    set(gca,"YMinorGrid","on")


    % Subplot 4: Qalpha vs. Depth
    subplot(1, 4, 4); % 2x2 grid, 4th subplot

    if ~isempty(Qalpha)
        for i = 1:length(Qalpha_depth_min)
            plot([Qalpha(i), Qalpha(i)], [Qalpha_depth_min(i), Qalpha_depth_max(i)], 'LineWidth', 3,'Marker','o'); hold on;
        end
        plot(0.8 * Qi,depth_index,'LineStyle','--','Color',[0.5, 0.5, 0.5],'LineWidth',3)
        plot(1.2 * Qi,depth_index,'LineStyle','--','Color',[0.5, 0.5, 0.5],'LineWidth',3)
        ylabel('Depth (km)');
        xlabel('Qalpha');
        title('Qalpha');
        grid on;
        set(gca,'YDir','reverse','FontSize',16,'FontWeight','bold')
        ylim([-10 1600])
        set(gca,"XMinorTick","on")
        set(gca,"YMinorTick","on")
        set(gca,"XMinorGrid","on")
    set(gca,"YMinorGrid","on")
    end
    
end
