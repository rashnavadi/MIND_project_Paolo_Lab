%% Section0: This script is going to call on event_timing_merge to merge event files per run and
% adjust the timing to be "absolute time". This means that the events from
% run 2 do not start at time 0, but at the end of run 1.

% This script will generate multiple text
% files with increasing numbers of events.

% Requirements:
% Main subject folders (named MinDch_Sub[number] for each subject)
%   folder called Event_timings
%       Event file descriptor as a .txt file  
%       Event timing files for each run, saved as Sub [subject number][EV] 
%       Run [run].txt (e.g., 'MIND_fMRI_Sub_034_EV1_Run1.txt')

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
    mkdir(fullfile(Paths.dataPath, SubjectFolders(CurrSub), '4_Analyses', 'Step2'))
    Subject.ID         = SubjectFolders(CurrSub);
    Paths.analysisPath = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step2');

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

        Paths.envelopesPath = fullfile(Paths.indexEventPath, 'Envelopes');
        if ~isfolder(Paths.envelopesPath)
            mkdir(Paths.envelopesPath)
        end
        
        Paths.binnedTimingsPath = fullfile(Paths.indexEventPath, 'BinnedTimings');
        if ~isfolder(Paths.binnedTimingsPath)
            mkdir(Paths.binnedTimingsPath)
        end
        
        clear tmp_*
    
        %% Section2: Calculate Envelope
    
        for run_looper = 1:numel(Subject.EEGruns)
            Run = Subject.EEGruns(run_looper);
            Run = extractAfter(Run,'Run');
            
            % Load Data 
            Internal.data_f_Path = fullfile(Paths.dataPath, Subject.ID, '2_EEG', '2_Cleaned', append(Subject.ID, '_Run', num2str(Run), '_Cleaned'));
            folder_dir = dir(Internal.data_f_Path);
            folder_dir = extractfield(folder_dir, 'name');
            folder_dir = string(folder_dir);
            folder_dir = folder_dir';
            i_fileName = contains(folder_dir, '.bin.');
            i_fileName = i_fileName * -1 + 1;
            i_fileName = find(i_fileName);
            Internal.data_f_Name = folder_dir(i_fileName);
            i_fileRemove = find(contains(Internal.data_f_Name, '.bin'));
            Internal.data_f_Name = folder_dir(i_fileRemove);
            Internal.log_f_Path = fullfile(Paths.dataPath, Subject.ID, '2_EEG', '2_Cleaned', append(Subject.ID, '_Run', num2str(Run),'_Cleaned'));
            Internal.log_f_Name = strcat(Subject.ID, '_Run', num2str(Run),'_ProcessingLog_&_ChannelLabels.txt');
    
            % Open Log File
            tmp_fid = fopen(fullfile(Internal.log_f_Path,Internal.log_f_Name));
            tmp_headerCell = textscan(tmp_fid,'%s', 'Delimiter', ':');
            tmp_headerCell = tmp_headerCell{1};
            fclose(tmp_fid);
    
            % Convert cell array to string array
            for tmp_looper = 1:length(tmp_headerCell)
                tmp_header(tmp_looper) = string([tmp_headerCell{tmp_looper}]);    
            end

            % Get Relevant Variables 
            EEG.subjectNumber = str2double(tmp_header((find(strcmp(tmp_header, 'subjectNumber')) + 1)));
            EEG.runNumber = str2double(tmp_header((find(strcmp(tmp_header, 'runNumber'),1) + 1)));
            EEG.s_Rate = str2double(tmp_header((find(strcmp(tmp_header, 's_Rate'),1) + 1)));
            EEG.n_Samples = str2double(tmp_header((find(strcmp(tmp_header, 'n_Samples'),1) + 1)));
            EEG.n_Channels = str2double(tmp_header((find(strcmp(tmp_header, 'n_Channels'),1) + 1)));
            tmp_i_first_c_Label = (find(strcmp(tmp_header, 'c_Labels (Column)'),1) + 1);
            EEG.c_Labels = tmp_header(tmp_i_first_c_Label:end);

            % Select Event File type and load
            Internal.event_f_Path = fullfile(Paths.dataPath, SubjectFolders(CurrSub), '2_EEG', '3_Events', Subject.Events(EV_loop), 'EnvelopeTimings');
            Internal.event_f_Names = strcat(Subject.ID, '_', Subject.EEGruns(run_looper), '_', Subject.Events(EV_loop), '.txt');
    
            % Open Event File
            tmp_fid = fopen(fullfile(Internal.event_f_Path, Internal.event_f_Names));
            tmp_events = fscanf(tmp_fid,'%f %f %f');
            fclose(tmp_fid);
            tmp_events = reshape(tmp_events, [3, (length(tmp_events)/3)]);
            tmp_events(2:3,:) = [];
            EEG.Events = round((tmp_events * EEG.s_Rate), 0);

            % Read EEG Channel Data
            tmp_fid = fopen(fullfile(Internal.data_f_Path, Internal.data_f_Name));
            EEG.c_Data = fread(tmp_fid, [EEG.n_Channels, EEG.n_Samples], 'float');  
            % =========================
            % Added by Tahereh, April, 2024
            % resize EEG.c_Data_Bipolar to EEG.c_Data (or n_Samples)
            EEG.c_Data = EEG.c_Data(:, 1:EEG.n_Samples);
            % =========================
            fclose(tmp_fid);

            % Set Parameters 
            % Set Display Channels
            % First Pass (Init Vars) of stripping out alpha prefix of contact name
            EEG.i_electrodeLastContacts = [];
            tmp_loops = 1;
            tmp_Label = char(EEG.c_Labels(1));
            tmp_alpha = isstrprop(tmp_Label, 'alpha');
            tmp_i_Alpha = find(tmp_alpha);  
            tmp_alphaCompStem = tmp_Label(tmp_i_Alpha);

            % Loop all other contact names, strip alpha prefix, save array indices 
            % where alpha prefix changes ("electrode ends")
            for tmp_looper = 2:length(EEG.c_Labels)
                    tmp_Label = char(EEG.c_Labels(tmp_looper));
                    tmp_alpha = isstrprop(tmp_Label, 'alpha');
                    tmp_i_Alpha = find(tmp_alpha);   
                    tmp_alphaStem = tmp_Label(tmp_i_Alpha);
                        if ((strlength(tmp_alphaStem) ~= strlength(tmp_alphaCompStem))||(string(tmp_alphaStem) ~= string(tmp_alphaCompStem)))
                            EEG.i_electrodeLastContacts(tmp_loops) = tmp_looper - 1;
                            tmp_loops = tmp_loops + 1;  
                            tmp_alphaCompStem = tmp_alphaStem; 
                        end
            end

            % Add last contact of last electrode end, because loop ends before
            % saving it
            EEG.i_electrodeLastContacts(tmp_loops) = length(EEG.c_Labels);

            % Make bipolars 
            tmp_loops = 1;
            for tmp_looper = 2:EEG.n_Channels
                if ~ismember((tmp_looper - 1), (EEG.i_electrodeLastContacts))
                    EEG.c_Data_Bipolar((tmp_loops),:) = EEG.c_Data((tmp_looper - 1),:) - EEG.c_Data(tmp_looper,:);
                    tmp_loops = tmp_loops + 1;
                end
            end

            % Re-label to bipolar 
            tmp_loops = 1;
            for tmp_looper = 2:EEG.n_Channels
                if ~ismember((tmp_looper - 1), (EEG.i_electrodeLastContacts))
                    EEG.c_Labels_Bipolar(tmp_loops) = append(EEG.c_Labels(tmp_looper - 1),'-',EEG.c_Labels(tmp_looper));
                    tmp_loops = tmp_loops + 1;
                end
            end

            % Select Channel
            if exist('SelectionMade', 'var') == 0
                    tmp_prompt = 'Please select the channel:';
                    tmp_mode = 'single';
                    tmp_title = 'Channel Selection';
                    [SelectionIndex, SelectionMade] = listdlg('ListString', string(EEG.c_Labels_Bipolar), 'SelectionMode', tmp_mode, 'Name', tmp_title, 'PromptString', tmp_prompt);
            elseif exist('SelectionMade', 'var') == 1
        
            end

            EEG.c_2Eval = SelectionIndex;
            c_Evaled = EEG.c_Labels_Bipolar(EEG.c_2Eval);
            tmp_c_data2Eval = EEG.c_Data_Bipolar(EEG.c_2Eval,:);

            % Set time and sample window 
            EEG.t_Window_pre = 200;
            EEG.t_Window_post = 300;
            EEG.s_Window_pre = (EEG.t_Window_pre/1000)*EEG.s_Rate;
            EEG.s_Window_post = (EEG.t_Window_post/1000)*EEG.s_Rate;  

            % Epoch Data
            EEG.Epoch = [];
            for tmp_looper = 1:length(EEG.Events)
                %========================
                % Added by Tahereh, April 2024
                if (EEG.Events(tmp_looper) + EEG.s_Window_post) > EEG.n_Samples
                    EEG.Events(tmp_looper) = EEG.n_Samples - EEG.s_Window_post;
                end
                %========================
                EEG.Epoch(tmp_looper,:) = tmp_c_data2Eval((EEG.Events(tmp_looper) - EEG.s_Window_pre):(EEG.Events(tmp_looper) + EEG.s_Window_post));
            end
            all(run_looper).Events = EEG.Epoch;

            % Calculate Stats
            signal_meanOfSamples = mean(EEG.Epoch);
            signal_SDofSamples = std(EEG.Epoch);

            %========================
            % Added by Tahereh, April 2024
            % if there is only one IED occurrence then skip the
            % normalization as it makes a 0/0 or Nan==> error
            if sum(signal_meanOfSamples ~= signal_meanOfSamples(1)) == 0  % check if all values in signal_meanOfSamples are the same
                norm_signal_meanOfSamples = signal_meanOfSamples;
                norm_signal_SDOfSamples = signal_SDofSamples;
            else
            %========================
                norm_signal_meanOfSamples = (signal_meanOfSamples - min(signal_meanOfSamples))/(max(signal_meanOfSamples) - min(signal_meanOfSamples));
                norm_signal_SDOfSamples = (signal_SDofSamples - min(signal_SDofSamples))/(max(signal_SDofSamples) - min(signal_SDofSamples));
                signal_SDofAllEpochs = sum(norm_signal_SDOfSamples);
                signal_meanOfAllEpochs = sum(norm_signal_meanOfSamples);
                AUC_Approx = sum(abs(2.*(norm_signal_meanOfSamples - norm_signal_SDOfSamples)));
            end
            % Save Average Spike 
            tmp_sampleFigure = figure;
            tmp_sampleFigure.Units = 'normalized';
            tmp_sampleFigure.Position = [0 0 1 1];
            tmp_sampleFigure.Visible = 'off';
            tmp_FigureName = 'Average IED Spike';
            tmp_sampleFigure.Name = tmp_FigureName;

            % Plot Sample Data 
            % =========================
            % Added by Tahereh, April, 2024
            % Check if tmp_events is empty or contains NaNs (i.e., no IED event is recorded)
            if isempty(tmp_events) || any(isnan(tmp_events))
                tmp_events_flag = false; % Set a flag to indicate that tmp_events is zero
                tmp_events = 0; % Set tmp_events to zero
            else
                tmp_events_flag = true; % Set a flag to indicate that tmp_events is not zero
            end

            if tmp_events_flag
            % =========================
                plot(signal_meanOfSamples, 'LineWidth', 3, 'Color', 'k')
                hold on 
                plot((signal_meanOfSamples - signal_SDofSamples), 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
                plot((signal_meanOfSamples + signal_SDofSamples), 'LineWidth', 3, 'Color',[0.6350 0.0780 0.1840])
                hold off
                xticks([1, 250, 500, 750, 1000, 1250])
                xticklabels({'-200', '-100', '0', '100', '200', '300'})
                xlabel('Time from IED Peak (ms)', 'FontSize', 18, 'FontWeight', 'bold')
            % =========================
            % Added by Tahereh, April, 2024
            if isscalar(norm_signal_SDOfSamples) % checks if norm_signal_SDOfSamples is a scalar (i.e., a single value)
                xlim([1, 2]);
            else
                % =========================
                xlim([1, length(norm_signal_SDOfSamples)])
                ylabel('Amplitude (uV)', 'FontSize', 18, 'FontWeight', 'bold')
                title(append(Subject.EEGruns(run_looper), ' Average IED Waveform with Envelope'), 'FontSize', 24, 'FontWeight', 'bold')
                xline(500, '--k', 'LineWidth', 2);
                ax = gca;
                ax.FontSize = 14;
            end
            % =========================
            % Added by Tahereh, April, 2024
            else
                % Handle the case when tmp_events is zero
                % For example, set xlim to [0, 1] if tmp_events is zero
                xlim([0, 1]);
            end
            % =========================

            %  Line(EEG.epochCentre, 'LineWidth', 2, 'Color', 'k');
            tmp_figFileName = append(Subject.ID,'_Run', num2str(Run), '_', Subject.Events(EV_loop), '_Average_IED_Spike'); 
            saveas(tmp_sampleFigure, fullfile(Paths.envelopesPath, tmp_figFileName), 'png')

            % Save Output .txt
            tmp_log_f_Name = append(Subject.ID,'_Run', string(EEG.runNumber), '_', Subject.Events(EV_loop),'_LogFile.txt'); 
            tmp_fid = fopen(fullfile(Paths.envelopesPath, tmp_log_f_Name), 'w');
            fprintf(tmp_fid, ['%s %s %s All Runs Envelope Log File \n',...
            'Channel Selected: %s \n', ... 
            'Pre-Event Window: %.0f ms \n', ...
            'Post-Event Window: %.0f ms \n', ... 
            'Sum of Signal Mean of All Epochs: %.0f \n', ...
            'Sum of Signal SD of All Epochs: %.0f \n', ...
            'Approximated AUC (normalized): %.0f \n', ...
            ], Subject.ID, Subject.EEGruns(run_looper), Subject.Events(EV_loop), EEG.c_Labels_Bipolar(EEG.c_2Eval), EEG.t_Window_pre, EEG.t_Window_post, ... 
            signal_meanOfAllEpochs, signal_SDofAllEpochs, AUC_Approx);
            fclose(tmp_fid);
            
            EEG.runNumber = str2double(tmp_header((find(strcmp(tmp_header, 'runNumber'),1) + 1)));
            EEG.s_Rate = str2double(tmp_header((find(strcmp(tmp_header, 's_Rate'),1) + 1)));
            EEG.n_Samples = str2double(tmp_header((find(strcmp(tmp_header, 'n_Samples'),1) + 1)));
   
            clear tmp_* EEG Internal
        end
    
        % Calculate full envelope stats
        combined_Events = all(1).Events;
        if numel(Subject.EEGruns) > 1
            for loops = 2:numel(Subject.EEGruns)
                combined_Events = vertcat(combined_Events, all(loops).Events);
            end
        end

        signal_meanOfSamples = mean(combined_Events);
        signal_SDofSamples = std(combined_Events);
        norm_signal_meanOfSamples = (signal_meanOfSamples - min(signal_meanOfSamples))/(max(signal_meanOfSamples) - min(signal_meanOfSamples));
        norm_signal_SDOfSamples = (signal_SDofSamples - min(signal_SDofSamples))/(max(signal_SDofSamples) - min(signal_SDofSamples));
        x = 1:length(norm_signal_SDOfSamples);
        signal_SDofAllEpochs = sum(norm_signal_SDOfSamples);
        signal_meanOfAllEpochs = sum(norm_signal_meanOfSamples);
        AUC_Approx = sum(abs(2.*(norm_signal_meanOfSamples - norm_signal_SDOfSamples)));
    
        % Save Average Spike 
        tmp_sampleFigure = figure;
        tmp_sampleFigure.Units = 'normalized';
        tmp_sampleFigure.Position = [0 0 1 1];
        tmp_sampleFigure.Visible = 'off';
        tmp_FigureName = 'Average IED Spike';
        tmp_sampleFigure.Name = tmp_FigureName;
    
        % Plot Sample Data 
        plot(x,signal_meanOfSamples, 'LineWidth', 3, 'Color', 'k')
        hold on 
        plot(x,(signal_meanOfSamples - signal_SDofSamples), 'LineWidth', 3, 'Color', [0 0 0.74])
        plot(x,(signal_meanOfSamples + signal_SDofSamples), 'LineWidth', 3, 'Color',[0.7 0 0])
        hold off
        xticks([1, 250, 500, 750, 1000, 1250])
        xticklabels({'-200', '-100', '0', '100', '200', '300'})
        xlabel('Time from IED Peak (ms)', 'FontSize', 18, 'FontWeight', 'bold')
        xlim([1, length(norm_signal_SDOfSamples)])
        ylabel('Amplitude (uV)', 'FontSize', 18, 'FontWeight', 'bold')
        title('Average IED Waveform with Envelope', 'FontSize', 24, 'FontWeight', 'bold')
        xline(500, '--k', 'LineWidth', 2);
        ax = gca;
        ax.FontSize = 14;
        tmp_figFileName = append(Subject.ID,'_AllRuns_Average_', Subject.Events(EV_loop), '_Spike'); 
        saveas(tmp_sampleFigure, fullfile(Paths.envelopesPath, tmp_figFileName), 'png')
    
        % Save Average Spike 
        tmp_sampleFigure = figure;
        tmp_sampleFigure.Units = 'normalized';
        tmp_sampleFigure.Position = [0 0 1 1];
        tmp_sampleFigure.Visible = 'off';
        tmp_FigureName = 'Average IED Spike';
        tmp_sampleFigure.Name = tmp_FigureName;
    
        % Plot Sample Data 
        patch([x, fliplr(x)], [(signal_meanOfSamples + signal_SDofSamples), fliplr((signal_meanOfSamples - signal_SDofSamples))], [0.71 0.91 1])
        hold on 
        plot(x,signal_meanOfSamples, 'LineWidth', 3, 'Color', 'k')
        plot(x,(signal_meanOfSamples - signal_SDofSamples), 'LineWidth', 3, 'Color', [0 0 0.74])
        plot(x,(signal_meanOfSamples + signal_SDofSamples), 'LineWidth', 3, 'Color',[0.7 0 0])
        hold off
        xticks([1, 250, 500, 750, 1000, 1250])
        xticklabels({'-200', '-100', '0', '100', '200', '300'})
        xlabel('Time from IED Peak (ms)', 'FontSize', 18, 'FontWeight', 'bold')
        xlim([1, length(norm_signal_SDOfSamples)])
        ylabel('Amplitude (uV)', 'FontSize', 18, 'FontWeight', 'bold')
        title('Average IED Waveform with Envelope', 'FontSize', 24, 'FontWeight', 'bold')
        xline(500, '--k', 'LineWidth', 2);
        ax = gca;
        ax.FontSize = 14;
        tmp_figFileName = append(Subject.ID,'_AllRuns_Average_', Subject.Events(EV_loop), '_Spike_Shaded'); 
        saveas(tmp_sampleFigure, fullfile(Paths.envelopesPath, tmp_figFileName), 'png')
        tmp_data_f_Name = fullfile(Paths.envelopesPath, append(Subject.ID,'_Envelope_Data.txt'));
        tmp_envelopeData(1,:) = (signal_meanOfSamples - signal_SDofSamples)';
        tmp_envelopeData(2,:) = (signal_meanOfSamples)';
        tmp_envelopeData(3,:) = (signal_meanOfSamples + signal_SDofSamples)';
        dlmwrite(tmp_data_f_Name, tmp_envelopeData);
        tmp_log_f_Name = append(Subject.ID,'_AllRuns_Average', Subject.Events(EV_loop), '_Envelope_LogFile.txt'); 
        tmp_fid = fopen(fullfile(Paths.envelopesPath, tmp_log_f_Name), 'w');
        fprintf(tmp_fid, ['%s %s Envelope Log File \n',...
        'Channel Selected: %s \n', ... 
        'Pre-Event Window: %.0f ms \n', ...
        'Post-Event Window: %.0f ms \n', ... 
        'Sum of Signal Mean of All Epochs: %.0f \n', ...
        'Sum of Signal SD of All Epochs: %.0f \n', ...
        'Approximated AUC (normalized): %.0f \n', ...
        ], Subject.ID, Subject.Events(EV_loop), c_Evaled, 200, 300, ... 
        signal_meanOfAllEpochs, signal_SDofAllEpochs, AUC_Approx);
        fclose(tmp_fid);
    
        %% Section3: read in the event timings for the appropriate runs.
        
        % Merge Events (Loop through Events)
        
        %---------------------------------------------------------------------------
        % Merge events from all runs for one subject and one EV at a time. 
        % --------------------------------------------------------------------------
        
        % Each new run is assigned a new line in the struct, even the empty ones, which will be
        % =========================
        % Added by Tahereh, April, 2024
        % this line added for the cases with more than one type of IED,
        % without the line below, the directory of Path.eventsPath is not
        % correct 
        Paths.eventsPath   = fullfile(Paths.dataPath, SubjectFolders(CurrSub), '2_EEG', '3_Events');
        % =========================
        Paths.eventsPath = fullfile(Paths.eventsPath, Subject.Events(EV_loop));
        eventFiles = dir(Paths.eventsPath);
        eventFiles = string(extractfield(eventFiles, 'name'));
        eventFiles = eventFiles';
        eventFiles_remove = contains(eventFiles, '._');
        eventFiles_remove = eventFiles_remove * -1 + 1;
        eventFiles_remove = find(eventFiles_remove);
        eventFiles = eventFiles(eventFiles_remove);
        eventFiles_select = find(contains(eventFiles, '.txt'));
        eventFiles = eventFiles(eventFiles_select);
        
        clear eventFiles_remove eventFiles_select
        
        % Read in event file events
        for i = 1:numel(eventFiles)
            % =====================
            % ADDED BY TAHEREH, APRIL 2024
            % Check if the file is empty
            if exist(fullfile(Paths.eventsPath, eventFiles(i)), 'file') == 0 || ... % check if the file exists
                    (exist(fullfile(Paths.eventsPath, eventFiles(i)), 'file') == 2 && ... % check if the file exist and is empty
                    numel(importdata(fullfile(Paths.eventsPath, eventFiles(i)))) == 0)
                continue; % Skip to the next iteration if the file is empty
            end

            % Read the file if it's not empty
            % =====================
            Event(i).timings = dlmread(fullfile(Paths.eventsPath, eventFiles(i)));
        end
        
        % Copy for timing manipulation
        for i = 1:numel(eventFiles)
            Event(i).absoluteTime = Event(i).timings;
        end
        
        % Get the number of volume per run from the design file
        for i = 1:numel(eventFiles)
            designFile       = strcat(Subject.fMRIruns(i), '.feat', '/design.fsf');
            Paths.designPath = fullfile(Paths.fMRIrunsPath, designFile);
            FSL_Design       = fopen(Paths.designPath);
            tmp_Design       = textscan(FSL_Design,'%s','delimiter','\n','whitespace','');
            fclose(FSL_Design);
            
            tmp_Design     = tmp_Design{1,1};
            n_volumes_line = find(contains(tmp_Design, 'set fmri(npts)'));
            n_volumes      = tmp_Design{n_volumes_line, 1};
            n_volumes      = n_volumes(16:end);
            n_volumes      = str2num(n_volumes);
            runTime        = n_volumes * 1.5;
            MergeEEG.runTime(i) = runTime;
        end
        
        MergeEEG.runTime = MergeEEG.runTime';
        
        % Change the zero time for all runs after run1 relative to time zero in run1
        % DOESN'T WORK FOR fMRI CORRECTED TIMINGS. 
        for j = 2:numel(eventFiles)
            if j == 2
                
                Added_time = MergeEEG.runTime(j-1);
                Event(j).absoluteTime(:,1) = (Event(j).timings(:,1) + Added_time);
                
            else
                Added_time = sum([MergeEEG.runTime(1:j-1)]);
                Event(j).absoluteTime(:,1) = (Event(j).timings(:,1) + Added_time);
            end
        end
        
        % Should there be a catch for empty files? hm. 
        
        % This part is going to make a new variable that is the vertical
        % combination of the timings

        % =====================
        % ADDED BY TAHEREH, APRIL 2024
          if numel(Event) == 1 % there is only one fmri run and one IED file
            Merged_event_timings = vertcat(Event(1).timings, Event(1).absoluteTime);
          end
        % =====================

        for k = 2:numel(Event)
          if k == 2
            Merged_event_timings = vertcat(Event(1).timings, Event(k).absoluteTime);
          else
            Merged_event_timings = vertcat(Merged_event_timings, Event(k).absoluteTime);
              
          end
        end
        
        % Finally, the merged event timings variable will be stored as a text file
        % and the total number of events will be saved as a seperate text tile
        
        output_filename = fullfile(Paths.binnedTimingsPath, append(Subject.ID, '_', Subject.Events(EV_loop), '_Merged_Timings.txt'));
        % =====================
        % ADDED BY TAHEREH, APRIL 2024
        % check the length of merged fMRI data, its duration might be
        % less than the time that IED has occurred so it gives error
        fMRI_data_filename = fullfile(Paths.dataPath, Subject.ID, '4_Analyses/Step1/MergedFunctionals/', strcat(Subject.ID, '_Merged_functional_data_', Subject.fMRIruns(1), '-', Subject.fMRIruns(end), '.nii.gz'));
        if exist(fMRI_data_filename)
            fMRI_data_filename = fMRI_data_filename;
        else % there is only one fMRI run and there is no merged functional file
        fMRI_data_filename = fullfile(Paths.dataPath, Subject.ID, '3_MRI/2_Functionals/2_Preprocessed', strcat(Subject.fMRIruns(1), '.feat'), strcat(Subject.ID,'_', Subject.fMRIruns(1),'_', 'filtered_func_ICA.nii.gz')); 
        end

        % define the fslinfo command
        fslinfoCommand = ['fslinfo ', char(fMRI_data_filename)];
        [status, output_fslinfo] = system(fslinfoCommand);
        dim4 ='';
        % Use regular expression to find and extract dim4 value
        expr = 'dim4\s+(\d+)'; % Regular expression pattern to find dim4
        matches = regexp(output_fslinfo, expr, 'tokens'); % Find matches
        dim4_value = str2double(matches{1}{1}); % Extract and convert to number

        events_to_keep = Merged_event_timings(:, 1) <= dim4_value;
        Merged_event_timings = Merged_event_timings(events_to_keep, :);
        disp(['The number of events that are occurring within the length of fMRI : ', num2str(size(Merged_event_timings))]);
        % =====================
        dlmwrite(output_filename, Merged_event_timings, 'delimiter', '\t');
        
        total_num_events = size(Merged_event_timings, 1);
        
        if total_num_events <150
            num_events_per_iteration = (10:10:total_num_events)';
            if ~isequal(num_events_per_iteration(end), total_num_events)
                num_events_per_iteration = vertcat(num_events_per_iteration, total_num_events);
            end
        elseif total_num_events >150 && total_num_events < 500
            num_events_per_iteration1 = (10:10:150)';
            num_events_per_iteration2 = (200:50:total_num_events)';
            num_events_per_iteration = vertcat(num_events_per_iteration1, num_events_per_iteration2);
            if ~isequal(num_events_per_iteration(end), total_num_events)
                num_events_per_iteration = vertcat(num_events_per_iteration, total_num_events);
            end
        elseif total_num_events >500 
            num_events_per_iteration1 = (10:10:150)';
            num_events_per_iteration2 = (200:50:500)';
            num_events_per_iteration3 = (600:100:total_num_events)';
            num_events_per_iteration = vertcat(num_events_per_iteration1, num_events_per_iteration2, num_events_per_iteration3);
            if ~isequal(num_events_per_iteration(end), total_num_events)
                num_events_per_iteration = vertcat(num_events_per_iteration, total_num_events);
            end
        end
        
        Num_events_per_iteration_name = fullfile(Paths.binnedTimingsPath, append(Subject.ID, '_', Subject.Events(EV_loop), '_Num_Events_Per_Iteration.txt'));
        dlmwrite(Num_events_per_iteration_name, num_events_per_iteration);
        
        % Once the number of events for each file are defined, we can then make a
        % variable with the appropriate number of event timings, then save these
        % timings to txt.
        for i = 1:numel(num_events_per_iteration)
            split_event_timing = Merged_event_timings(1:num_events_per_iteration(i),:,:);
            num_events_str = num2str(num_events_per_iteration(i));
            output_filename = fullfile(Paths.binnedTimingsPath, append(Subject.ID, '_', Subject.Events(EV_loop), '_', num_events_str, '_events.txt'));
            dlmwrite(output_filename, split_event_timing, 'delimiter', '\t');
            
            %This is going to keep track of the timing of the last event in each
            %file.
            tmp_last_event_timing_row = split_event_timing(end,:,:);
            Last_event_timing_log(i) = tmp_last_event_timing_row(1);
            
        end
        
        % This writes the last event timings to text for use later.
        Last_event_timing_log = Last_event_timing_log';
        Last_event_timing_log_filename = fullfile(Paths.binnedTimingsPath, append(Subject.ID, '_', Subject.Events(EV_loop), '_Last_Event_Timing.txt'));
        dlmwrite(Last_event_timing_log_filename, Last_event_timing_log);
        clear SelectionMade;
        % =====================
        % ADDED BY TAHEREH, APRIL 2024
        clear Last_event_timing_log % you need to clear this file otherwise it combines timings of the previous subject with the current one
    end % End event loop
    disp(append('Done ', Subject.ID));
    clear SelectionMade;
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