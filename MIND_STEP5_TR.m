%% Section0: ***** This script runs fsl, but only if matlab is opened through terminal using 
% /Applications/Matlab_R2019a.app/bin/matlab

%%% THIS SCRIPT takes ALPHASIM output and compares the result to the
%%% analysis including all events.

% REQUIREMENTS
%Folder - subject specific (e.g., MinDch_Sub034)
%   Folder - Merged_runs
%       Folder - Feat results for each number of events
%           Folder - stats also containing alphasim results
%   Folder - Event_timings
%       Files - total number of events per iteration - for naming purposes.

%% Section1: Create batch for execution

% Create initial project paths
Paths.script             = mfilename( 'fullpath');
[Paths.scriptPath, ~, ~] = fileparts(Paths.script);
Paths.dataPath           = extractBefore(Paths.scriptPath,"/4_Scripts");

% Warn the user that Matlab must be initiated through terminal in order for
% the script to work
tmp_Box.WindowStyle = 'non-modal';
tmp_Box.Interpreter = 'tex';
waitfor(msgbox({'\fontsize{16} MIND must be run in Matlab initialized in terminal.'; ...
    '\fontsize{16} (/Applications/MATLAB\_R20##.app/bin/matlab &)'}, 'MIND', 'help', tmp_Box));

% Find all potential subjects for analysis
SubjectFolders = dir(Paths.dataPath);
SubjectFolders = string(extractfield(SubjectFolders, 'name'));
SubjectFolders = SubjectFolders';
tmp_SubjectFolders_Relevant = find(contains(SubjectFolders, 'MIND0', 'IgnoreCase', true));
SubjectFolders = SubjectFolders(tmp_SubjectFolders_Relevant);
tmp_Remove = contains(SubjectFolders, '._');
tmp_keep = find(tmp_Remove * -1 + 1);
SubjectFolders = SubjectFolders(tmp_keep);

% Generate GUI for the user to select the subjects they want to run
tmp_GUI_spacing = (30*numel(SubjectFolders)):-30:25;
tmp_GUI_height  = tmp_GUI_spacing(1) + 55;
tmp_GUI_spacing = tmp_GUI_spacing + 20;

tmp_fig = uifigure;
tmp_fig.Visible = 'on';
set(0,'units','pixels')
tmp_Pix_SS = get(0,'screensize');
tmp_figCtr = tmp_fig.Position(3:4)./2;
tmp_fig.Position(4) = 600;
tmp_fig.Position(3) = 350;
tmp_fig.Position(1:2) = (tmp_Pix_SS(3:4)./2) - tmp_figCtr;
tmp_fig.Name = 'MIND Subject Selection';
tmp_fig.Scrollable = 'on';
tmp_fig.Resize = 'on';

for tmp_loop = 1:numel(SubjectFolders)
   tmp_sub(tmp_loop) = uicheckbox(tmp_fig);
end
for tmp_loop = 1:numel(SubjectFolders)  
   tmp_sub(tmp_loop).Text = SubjectFolders(tmp_loop);
end
for tmp_loop = 1:numel(SubjectFolders)   
   tmp_sub(tmp_loop).Position(4) = 15;
end
for tmp_loop = 1:numel(SubjectFolders)  
   tmp_sub(tmp_loop).Position(3) = 75;
end
for tmp_loop = 1:numel(SubjectFolders)   
   tmp_sub(tmp_loop).Position(1) = 75;
end
for tmp_loop = 1:numel(SubjectFolders)   
   tmp_sub(tmp_loop).Position(2) = tmp_GUI_spacing(tmp_loop);
end
for tmp_loop = 1:numel(SubjectFolders)   
   tmp_sub(tmp_loop).Value = 0;
end

tmp_btn = uibutton(tmp_fig, 'ButtonPushedFcn', {@AddSubject, tmp_sub});
tmp_btn.Text = 'DONE';
tmp_btn.Position = [175 (tmp_GUI_height/2) 100 25];

tmp_fig.Visible = 'on';

waitfor(tmp_btn,'Enable', 'off')
close(tmp_fig)

tmp_SubjectStates = tmp_SubjectStates';

% Remove subjects not selected from the GUI
tmp_SubjectStates  = string(tmp_SubjectStates);
tmp_RemoveSubjects = find(contains(tmp_SubjectStates, 'true', 'IgnoreCase', true));
SubjectFolders     = SubjectFolders(tmp_RemoveSubjects);

clear tmp_*

for CurrSub = 1:numel(SubjectFolders)
    % Initialize and set subject variables
    mkdir(fullfile(Paths.dataPath, SubjectFolders(CurrSub), '4_Analyses', 'Step5'))
    Subject.ID         = SubjectFolders(CurrSub);
    Paths.analysisPath = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step5');

    Paths.fMRIrunsPath     = fullfile(Paths.dataPath, SubjectFolders(CurrSub), '3_MRI', '2_Functionals', '2_Preprocessed');
    tmp_fMRIrunPath        = dir(Paths.fMRIrunsPath);
    tmp_fMRIrunPath        = string(extractfield(tmp_fMRIrunPath, 'name'));
    tmp_fMRIrunPath        = tmp_fMRIrunPath';
    tmp_fMRIrunPath_Remove = find(contains(tmp_fMRIrunPath, 'Run'));
    tmp_fMRIrunPath        = tmp_fMRIrunPath(tmp_fMRIrunPath_Remove);
    tmp_fMRIrunPath_Remove = contains(tmp_fMRIrunPath, 'Exclude');
    tmp_fMRIrunPath_Remove = find(tmp_fMRIrunPath_Remove * -1 + 1);
    tmp_fMRIrunPath        = tmp_fMRIrunPath(tmp_fMRIrunPath_Remove);
    tmp_fMRIrunPath        = extractBefore(tmp_fMRIrunPath, '.feat');
    Subject.fMRIruns       = tmp_fMRIrunPath;

    Paths.EEGrunsPath     = fullfile(Paths.dataPath, SubjectFolders(CurrSub), '2_EEG', '2_Cleaned');
    tmp_EEGrunPath        = dir(Paths.EEGrunsPath);
    tmp_EEGrunPath        = string(extractfield(tmp_EEGrunPath, 'name'));
    tmp_EEGrunPath        = tmp_EEGrunPath';
    tmp_EEGrunPath_Remove = find(contains(tmp_EEGrunPath, 'Run'));
    tmp_EEGrunPath        = tmp_EEGrunPath(tmp_EEGrunPath_Remove);
    tmp_EEGrunPath_Remove = contains(tmp_EEGrunPath, 'Exclude');
    tmp_EEGrunPath_Remove = find(tmp_EEGrunPath_Remove * -1 + 1);
    tmp_EEGrunPath        = tmp_EEGrunPath(tmp_EEGrunPath_Remove);
    tmp_EEGrunPath        = extractBefore(tmp_EEGrunPath, '_Cleaned');
    tmp_EEGrunPath        = extractAfter(tmp_EEGrunPath,'_');
    Subject.EEGruns       = tmp_EEGrunPath;

    Paths.eventsPath   = fullfile(Paths.dataPath, SubjectFolders(CurrSub), '2_EEG', '3_Events');
    tmp_EventPath        = dir(Paths.eventsPath);
    tmp_EventPath        = string(extractfield(tmp_EventPath, 'name'));
    tmp_EventPath        = tmp_EventPath';
    tmp_EventPath_Remove = find(contains(tmp_EventPath, 'IED'));
    tmp_EventPath        = tmp_EventPath(tmp_EventPath_Remove);
    tmp_EventPath_Remove = contains(tmp_EventPath, 'MIND');
    tmp_EventPath_Remove = find(tmp_EventPath_Remove * -1 + 1);
    tmp_EventPath        = tmp_EventPath(tmp_EventPath_Remove);
    Subject.Events       = tmp_EventPath;

    for EV_loop = 1:numel(Subject.Events)

        % Make output folders
        Paths.indexEventPath = fullfile(Paths.analysisPath, Subject.Events(EV_loop));
        if ~isfolder(Paths.indexEventPath)
            mkdir(Paths.indexEventPath)
        end

        Paths.screenshotsPath = fullfile(Paths.indexEventPath, 'MapScreenshots');
        if ~isfolder(Paths.screenshotsPath)
            mkdir(Paths.screenshotsPath)
        end
        
        clear tmp_*
        
        %% Section2: Run Stats
        
        % This part finds the number of events per iteration, for naming purposes
        num_events_per_iteration_filename = fullfile(Paths.dataPath, Subject.ID, '4_Analyses/Step2', Subject.Events(EV_loop), 'BinnedTimings', strcat(Subject.ID, '_', Subject.Events(EV_loop), '_Num_Events_Per_Iteration.txt'));
        num_events_per_iteration = dlmread(num_events_per_iteration_filename);
        
        % This part creates a binarized mask of voxels for the total number of events
        
        amag_folder = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step4', Subject.Events(EV_loop), strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(end)), '_events'), 'Amalgamated');
        target_Clusters_path = fullfile(amag_folder, append(Subject.Events(EV_loop), '_', string(num_events_per_iteration(end)), '_events', '_zstat1_Composite.nii.gz')); 
        target_Clusters = niftiread(target_Clusters_path);
        
        % Binarize cluster map for all events 
        for tmp_x_Looper = 1:size(target_Clusters, 1)
            for tmp_y_Looper = 1:size(target_Clusters, 2)
                for tmp_z_Looper = 1:size(target_Clusters, 3)
                    if target_Clusters(tmp_x_Looper, tmp_y_Looper, tmp_z_Looper) > 0
                        target_Clusters(tmp_x_Looper, tmp_y_Looper, tmp_z_Looper) = 2; % Setting to 2 to allow true/false positive/negative differentiation when added to inerative maps.
                    end
                end
            end           
        end
        
        % This part finds the anatomical file for images for screenshot
        anat_filename = fullfile(Paths.dataPath, Subject.ID, '3_MRI', '1_Anatomicals', '2D_Anat', append(Subject.ID, '_2Danat_brain.nii.gz'));
        
        % This part is going to loop through all iterations of number of events 
        % and calculate contingency.
        
        for event_iteration = 1:(size(num_events_per_iteration, 1)-1)
            outfile = fullfile(Paths.indexEventPath, 'MapScreenshots', append(Subject.ID, "_", Subject.Events(EV_loop), '_', string(num_events_per_iteration(event_iteration)),...
                '_events.png'));
            
            iteration_cluster_map_path_Screenshot = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step4', Subject.Events(EV_loop), strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(event_iteration)), '_events'), 'Amalgamated', ...
                append(Subject.Events(EV_loop), '_', string(num_events_per_iteration(event_iteration)), '_events', '_zstat1_Composite_thresh_to2Danat.nii.gz'));

            iteration_cluster_map_path = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step4', Subject.Events(EV_loop), strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(event_iteration)), '_events'), 'Amalgamated', ...
                append(Subject.Events(EV_loop), '_', string(num_events_per_iteration(event_iteration)), '_events', '_zstat1_Composite_thresh.nii.gz'));

            %%% INSERT SAVE SCREENSHOT

            % Ensure the outfile has the correct extension
            outfile = strrep(outfile, '.nii.gz', '.png');  % Change the extension to .png
            % for lower versions of fsleyes uncomment this line and comment out
            % the one after
            %         screenshot_command = append("fsleyes render --outfile ", outfile, " -s lightbox ", anat_filename, " ",...
            %             iteration_cluster_map_path_Screenshot, " -zx 3 --sliceSpacing 8 --zrange 5 100 --ncols 5 --nrows 3 -cm red-yellow");

            % for fsleyes versioin like FSLeyes version 1.12.4 and probably
            % higher versions the line below works
            screenshot_command = append("fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 0.1 --zrange 0 1 --ncols 5 --nrows 3 --outfile ", outfile, " ", anat_filename, " ", iteration_cluster_map_path_Screenshot, " -cm red-yellow");
            system(screenshot_command);
            %%%% end save screenshot
            
            iteration_Clusters = niftiread(iteration_cluster_map_path);
            % Binarize cluster map for given iteration of events 
            for tmp_x_Looper = 1:size(iteration_Clusters, 1)
                for tmp_y_Looper = 1:size(iteration_Clusters, 2)
                    for tmp_z_Looper = 1:size(iteration_Clusters, 3)
                        if iteration_Clusters(tmp_x_Looper, tmp_y_Looper, tmp_z_Looper) > 0
                            iteration_Clusters(tmp_x_Looper, tmp_y_Looper, tmp_z_Looper) = 1; % Setting to 1 to allow true/false positive/negative differentiation when added to inerative maps.
                        end
                    end
                end           
            end
            % Make Confusion matrix (True/False Postive/Negative)
            confusion_values = iteration_Clusters + target_Clusters;
            % Split Confusion matrix (True/False Postive/Negative)
            n_true_positive = 0;
            n_false_positive = 0;
            n_true_negative = 0;
            n_false_negative = 0;
            for tmp_x_Looper = 1:size(confusion_values, 1)
                for tmp_y_Looper = 1:size(confusion_values, 2)
                    for tmp_z_Looper = 1:size(confusion_values, 3)
                        switch confusion_values(tmp_x_Looper, tmp_y_Looper, tmp_z_Looper)
                            case 3
                                n_true_positive = n_true_positive +1;
                            case 2
                                n_false_negative = n_false_negative + 1;
                            case 1
                                n_false_positive = n_false_positive + 1;
                            case 0
                                n_true_negative = n_true_negative + 1;
                        end
                    end
                end           
            end
            % GET CONTINGENCY STATISTICS: F1 SCORE & ASSOCIATED METRICS
            out_Table.n_Events(event_iteration) = num_events_per_iteration(event_iteration);
            out_Table.n_true_positive(event_iteration) = n_true_positive;
            out_Table.n_false_negative(event_iteration) = n_false_negative;
            out_Table.n_false_positive(event_iteration) = n_false_positive;
            out_Table.n_true_negative(event_iteration) = n_true_negative;
            out_Table.Accuracy(event_iteration) = (n_true_positive + n_true_negative)/(n_true_positive + n_true_negative + n_false_positive + n_false_negative);
            out_Table.Sensitivity(event_iteration) = n_true_positive/(n_true_positive + n_false_negative);
            out_Table.Specificity(event_iteration) = n_true_negative/(n_true_negative + n_false_positive);
            out_Table.Percision(event_iteration) = n_true_positive/(n_true_positive + n_false_positive);
            out_Table.F1(event_iteration) = (2 * (out_Table.Percision(event_iteration) * out_Table.Sensitivity(event_iteration))/(out_Table.Percision(event_iteration) + out_Table.Sensitivity(event_iteration)));
            out_Table.F05(event_iteration) = (1.25 * (out_Table.Percision(event_iteration) * out_Table.Sensitivity(event_iteration))/(0.25 * out_Table.Percision(event_iteration) + out_Table.Sensitivity(event_iteration)));
            out_Table.F2(event_iteration) = (5 * (out_Table.Percision(event_iteration) * out_Table.Sensitivity(event_iteration))/(4 * out_Table.Percision(event_iteration) + out_Table.Sensitivity(event_iteration)));
        end

        %%% INSERT SAVE SCREENSHOT # 2 - last one
        event_iteration = size(num_events_per_iteration, 1);
        iteration_cluster_map_path_Screenshot = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step4', Subject.Events(EV_loop), strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(end)), '_events'), 'Amalgamated', ...
            append(Subject.Events(EV_loop), '_', string(num_events_per_iteration(event_iteration)), '_events', '_zstat1_Composite_thresh_to2Danat.nii.gz'));
        outfile = fullfile(Paths.indexEventPath, 'MapScreenshots', strcat(Subject.ID, '_', Subject.Events(EV_loop), '_', string(num_events_per_iteration(event_iteration)),...
            '_events.png'));

        % for lower versions of fsleyes uncomment this line and comment out
        % the one after
        %         screenshot_command = append("fsleyes render --outfile ", outfile, " -s lightbox ", anat_filename, " ",...
        %             iteration_cluster_map_path_Screenshot, " -zx 3 --sliceSpacing 8 --zrange 5 100 --ncols 5 --nrows 3 -cm red-yellow");

        % for fsleyes versioin like FSLeyes version 1.12.4 and probably
        % higher versions the line below works
        screenshot_command = append("fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 0.1 --zrange 0 1 --ncols 5 --nrows 3 --outfile ", outfile, " ", anat_filename, " ", iteration_cluster_map_path_Screenshot, " -cm red-yellow");
        system(screenshot_command);
        %%%% End save screenshot

        % Save STATISTICS
        numevents = num_events_per_iteration(1:end-1);
        out_Table_tmp = structfun(@transpose,out_Table,'UniformOutput',false);
        stats_table = struct2table(out_Table_tmp);
        stats_table_name = fullfile(Paths.indexEventPath, strcat(Subject.ID, '_', Subject.Events(EV_loop), '_stats_table.xlsx'));
        writetable(stats_table, stats_table_name);
    end % End event loop
    disp(append('Done ', Subject.ID));
end % End subject loop

clear

%% Define Minor Functions for GUI to work %%

function AddSubject(src,event, tmp_sub)
    for looper = 1:length(tmp_sub)
        states(looper) = tmp_sub(looper).Value;
    end
    assignin('base','tmp_SubjectStates', states);
    tmp_btn = evalin('base', 'tmp_btn');
    tmp_btn.Enable = 'off';
    assignin('base','tmp_btn', tmp_btn);
    return
end
