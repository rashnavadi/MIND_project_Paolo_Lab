%% Section0: ***** This script runs fsl, but only if matlab is opened through terminal using 
% /Applications/Matlab_R2019a.app/bin/matlab

%%% THIS SCRIPT WILL RUN make_fsf.m to create DESIGN FILES FOR EACH LENGTH OF FMRI. 

% Requirements
% Folder - HRFs need to be stored somewhere - fill in location below. 
% Folder - Subject specific (e.g., Min_Dch_034)
%   File - design_template1 (general design, placeholders included
%            in locations requiring input - e.g., # volumes, data names)
%   Folder - Merged_runs
%       Files - one file, trimmed for each number of events, for each EV
%       File - volumes_per_iteration_per_run from 'make fsf files' script
%   Folder - Event timings
%       File - num_events_per_iteration from 'event timings separate' script
%   Folder - design_files (this will be where all the design files end up)

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
    mkdir(fullfile(Paths.dataPath, SubjectFolders(CurrSub), '4_Analyses', 'Step4'))
    Subject.ID         = SubjectFolders(CurrSub);
    Paths.analysisPath = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step4');

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

    clear tmp_*

    for EV_loop = 1:numel(Subject.Events)

        % Make output folders
        Paths.indexEventPath = fullfile(Paths.analysisPath, Subject.Events(EV_loop));
        if ~isfolder(Paths.indexEventPath)
            mkdir(Paths.indexEventPath)
        end

        img_dms = [64, 64, 24]; % default value - only change if necessary
        
        % Path to flobs
        flobs = ["peak3s", "peak5s", "peak7s", "peak9s"];
        Paths.flobPath = string(1:numel(flobs))';
        for event_looper = 1 : numel(flobs)
            flobPath = fullfile(Paths.scriptPath, 'HRFs', strcat(flobs(event_looper), '.flobs'));
            Paths.flobPath(event_looper) = flobPath;
        end
        
        % Path to template 
        Paths.templatePath = fullfile(Paths.scriptPath, 'Templates');
        
        %% Section2: Make Design Files 
        
        % Make an fsf file for different lengths of fMRI for minumum discharges
        % specifically, this will change the number of volumes, , number of pixels
        % (X x Y x Z x volumes), the 4D data, EVs, and HRF used to convolve EVs
        
        % This part loads in both the number of events and the number of volumes
        %%for each split file
        
        num_events_per_iteration_filename = fullfile(Paths.dataPath, Subject.ID, '4_Analyses/step2', Subject.Events(EV_loop), 'BinnedTimings', strcat(Subject.ID, '_', Subject.Events(EV_loop), '_Num_Events_Per_Iteration.txt'));
        num_events_per_iteration = dlmread(num_events_per_iteration_filename);
        
        volumes_per_iteration_name = fullfile(Paths.dataPath, Subject.ID, '4_Analyses/Step3', Subject.Events(EV_loop), append(Subject.ID, '_', Subject.Events(EV_loop), '_volumes_per_iteration.txt'));
        volumes_per_iteration = dlmread(volumes_per_iteration_name); 
        
        % This part loads the design template
        
        for event_looper = 1:numel(num_events_per_iteration)
            
            % Make iterative events output folder
            Paths.iterativeEventsPath = fullfile(Paths.indexEventPath, strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(event_looper)), '_events'));
            if ~isfolder(Paths.iterativeEventsPath)
                mkdir(Paths.iterativeEventsPath)
            end
        
            for flob_looper = 1:numel(flobs)
        
                 % Make iterative flobs output folder
                Paths.iterativeFlobsPath = fullfile(Paths.indexEventPath, strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(event_looper)), '_events'), flobs(flob_looper));
                if ~isfolder(Paths.iterativeFlobsPath)
                    mkdir(Paths.iterativeFlobsPath)
                end
        
                % Will need this throughout loop iteration
                num_volumes = num2str(volumes_per_iteration(event_looper));
                num_events = num2str(num_events_per_iteration(event_looper));
                
                % Open design template, change certain fields
                tmplt_ID = fopen(fullfile(Paths.templatePath,'MIND_Design_Template1.fsf'));
        
                FSL_tmplt = textscan(tmplt_ID,'%s','delimiter','\n','whitespace','');
        
                fclose(tmplt_ID);
                
                % Change number of volumes (changes for every iteration)
                FSL_tmplt{1,1} = regexprep(FSL_tmplt{1,1},'N_VOLUMES', num_volumes);
               
                % Change the total number of voxels in the analysis
                num_voxels = num2str(img_dms(1) * img_dms(2) * img_dms(3) * volumes_per_iteration(event_looper));
                FSL_tmplt{1,1} = regexprep(FSL_tmplt{1,1},'TOTAL_N_VOXELS', num_voxels);
               
                % Change Output folder 
                fMRI_outName = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step4', Subject.Events(EV_loop), append(Subject.Events(EV_loop), '_', num_events, '_events'), flobs(flob_looper));
                FSL_tmplt{1,1} = regexprep(FSL_tmplt{1,1},'OUT_FILE',fMRI_outName);
              
                % Change fmri filename 
                fMRI_filename = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step3', Subject.Events(EV_loop), append(Subject.ID, '_', Subject.Events(EV_loop), '_', num_events, '_events', '.nii.gz'));
                FSL_tmplt{1,1} = regexprep(FSL_tmplt{1,1},'fMRI_FILE',fMRI_filename);
             
                % Change EV file                    
                EV_filename = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'step2', Subject.Events(EV_loop), 'BinnedTimings', append(Subject.ID, '_', Subject.Events(EV_loop), '_', num_events, '_events.txt'));
                FSL_tmplt{1,1} = regexprep(FSL_tmplt{1,1},'EV_FILE', EV_filename);
               
                % Change HRF File 
                HRF_filename = fullfile(Paths.flobPath(flob_looper), 'hrfbasisfns.txt');
                FSL_tmplt{1,1} = regexprep(FSL_tmplt{1,1},'HRF_FILE', HRF_filename);
        
                % Make a new file, save it
                design_file_name = fullfile(Paths.iterativeFlobsPath, append(Subject.ID, '_', Subject.Events(EV_loop), '_design_', num_events, '_events', '_', flobs(flob_looper), '.fsf'));
                fileID = fopen(design_file_name,'wt', 'b', 'ISO-8859-1');
        
                for ii=1:size(FSL_tmplt{1,1},1)        
                    fprintf(fileID, '%s\n', FSL_tmplt{1,1}{ii,1});        
                end        
                fclose(fileID);
            end
        end
        
        %% Section3: Run fMRI analysis
        
        % This part loads in the number of events for naming purposes
        num_events_per_iteration_filename = fullfile(Paths.dataPath, Subject.ID, '4_Analyses/step2', Subject.Events(EV_loop), 'BinnedTimings', strcat(Subject.ID, '_', Subject.Events(EV_loop), '_Num_Events_Per_Iteration.txt'));
        num_events_per_iteration = dlmread(num_events_per_iteration_filename);
        
        % This part loops through each of the fsf's generated in 'make fsf files' 
        for i = 1:numel(num_events_per_iteration)
           for j = 1:numel(flobs)
                num_events = num2str(num_events_per_iteration(i));
                iterativeFlobsPath  = fullfile(Paths.indexEventPath, strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(i)), '_events'), flobs(j));
                iterativeDesignFile = append(Subject.ID, '_', Subject.Events(EV_loop), '_design_', num_events, '_events', '_', flobs(j), '.fsf');
                iterativePathFile   = fullfile(iterativeFlobsPath, iterativeDesignFile);
                feat_command = append('feat ', iterativePathFile);
                system(feat_command);
           end
        end
        
        %% Section4: Amalgamate Flobs
        
        for i = 1:numel(num_events_per_iteration)
            amag_folder = fullfile(Paths.indexEventPath, strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(i)), '_events'), 'Amalgamated');
        
            if ~isfolder(amag_folder)
                mkdir(amag_folder)
            end
        end
        % =====================
        % ADDED BY TAHEREH, APRIL 2024
        % Get Thresholded Images
        for looper = 1:numel(num_events_per_iteration)
            path = fullfile(Paths.indexEventPath, strcat(Subject.Events(EV_loop), '_', num2str(num_events_per_iteration(looper)), '_events'));
            % Flob 1 (peak3s)
            tmp_flob1_filepath = fullfile(path, 'peak3s.feat', 'thresh_zstat1.nii.gz');
            if exist(tmp_flob1_filepath)
                tmp_flob1 = niftiread(tmp_flob1_filepath);
                tmp_info1 = niftiinfo(tmp_flob1_filepath);
            else
                disp('threshold_zstate1.nii.gz for flob1 (3s) does not exist')
            end
            % Flob 2 (peak5s)
            tmp_flob2_filepath = fullfile(path, 'peak5s.feat', 'thresh_zstat1.nii.gz');
            if exist(tmp_flob2_filepath)
                tmp_flob2 = niftiread(tmp_flob2_filepath);
                tmp_info2 = niftiinfo(tmp_flob2_filepath);
            else
                disp('threshold_zstate1.nii.gz for flop 2 (5s) does not exist')
            end
            % Flob 3 (peak7s)
            tmp_flob3_filepath = fullfile(path, 'peak7s.feat', 'thresh_zstat1.nii.gz');
            if exist(tmp_flob3_filepath)
                tmp_flob3 = niftiread(tmp_flob3_filepath);
                tmp_info3 = niftiinfo(tmp_flob3_filepath);
            else
                disp('threshold_zstate1.nii.gz for flob3 (7s) does not exist')
            end
            % Flob 4 (peak9s)
            tmp_flob4_filepath = fullfile(path, 'peak9s.feat', 'thresh_zstat1.nii.gz');
            if exist(tmp_flob4_filepath)
                tmp_flob4 = niftiread(tmp_flob4_filepath);
                tmp_info4 = niftiinfo(tmp_flob4_filepath);
            else
                disp('threshold_zstate1.nii.gz for flob4 (9s) does not exist')
            end
            % Make New Nifti 
            tmp_out = zeros(size(tmp_flob1));
            % Get the max voxel 
            for tmp_lps = 1:size(tmp_flob1, 1)
                for tmp_ls = 1:size(tmp_flob1, 2)
                    for tmp_loop = 1:size(tmp_flob1, 3)
                        if exist ('tmp_flob1', 'var')
                            voxel(1) = tmp_flob1(tmp_lps, tmp_ls, tmp_loop);
                        else
                            voxel(1) = 0;
                        end
                        if exist ('tmp_flob2', 'var')
                            voxel(2) = tmp_flob2(tmp_lps, tmp_ls, tmp_loop);
                        else
                            voxel(2) = 0;
                        end
                        if exist('tmp_flob3', 'var')
                            voxel(3) = tmp_flob3(tmp_lps, tmp_ls, tmp_loop);
                        else
                            voxel(3) = 0;
                        end
                        if exist('tmp_flob4', 'var')
                            voxel(4) = tmp_flob4(tmp_lps, tmp_ls, tmp_loop);
                        else
                            voxel(4) = 0;
                        end
                        tmp_out(tmp_lps, tmp_ls, tmp_loop) = max(voxel);
                    end
                end           
            end
            tmp_out            = single(tmp_out);
            tmp_out_FileName   = append(Subject.Events(EV_loop), '_', string(num_events_per_iteration(looper)), '_events', '_zstat1_Composite.nii');
            tmp_out_FilePath   = fullfile(path, 'Amalgamated');
            tmp_out_PathFile   = fullfile(tmp_out_FilePath, tmp_out_FileName);
            tmp_info1.Filename = tmp_out_PathFile;
            niftiwrite(tmp_out, tmp_out_PathFile, tmp_info1);
            gzip(tmp_out_PathFile);
            delete(tmp_out_PathFile);
            cd(tmp_out_FilePath);
            Input1 = strcat(Subject.Events(EV_loop), '_', string(num_events_per_iteration(looper)), '_events', '_zstat1_Composite_ClusterIndex');
            Input2 = strcat(Subject.Events(EV_loop), '_', string(num_events_per_iteration(looper)), '_events', '_zstat1_Composite_thresh');
            Input3 = strcat(Subject.Events(EV_loop), '_', string(num_events_per_iteration(looper)), '_events', '_zstat1_Composite_Max.txt');
            Input4 = strcat(Subject.Events(EV_loop), '_', string(num_events_per_iteration(looper)), '_events', '_zstat1_Composite_ClusterTable.txt');
            % MODIFIED BY TAHEREH IN AUGUST 2024, I INCREASED THRESHOLD
            % FROM 3.1 TO 3.5 TO REDUCE NOISE AND HAVE MORE STATISTICALLY
            % SIGNIFICANT VALUES IN THE VOXELS
            tmp_command = append('cluster -i ', tmp_out_FileName, ' -t 3.5 --mm -o ', Input1, ' --othresh=', Input2, ' --olmax=', Input3, ' > ', Input4);
            system(tmp_command);
            % Execute the command
            system(tmp_command);
        % =====================
%             clear tmp_* 
        
            % Clean up design file
            for flob_looper = 1:numel(flobs)
                eventFolder  = strcat(Subject.Events(EV_loop), '_', string(num_events_per_iteration(looper)), '_events');
                flobFolderS  = flobs(flob_looper);
                flobFolderD  = strcat(flobs(flob_looper), '.feat');
                fileName     = strcat(Subject.ID, '_', Subject.Events(EV_loop), '_design_', string(num_events_per_iteration(looper)), '_events_', flobs(flob_looper), '.fsf');
                fileToMoveS  = fullfile(Paths.indexEventPath, eventFolder, flobFolderS, fileName);
                fileToMoveD  = fullfile(Paths.indexEventPath, eventFolder, flobFolderD);
            end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Register from func space to 2Danat and 3Danat space
            flirt_input    = fullfile(tmp_out_FilePath, append(Input2, '.nii.gz')); % composite image in func space
            flirt_2Danat   = fullfile(Paths.dataPath, Subject.ID, '3_MRI', '1_Anatomicals', '2D_Anat', strcat(Subject.ID, '_2Danat_brain.nii.gz'));
            flirt_3Danat   = fullfile(Paths.dataPath, Subject.ID, '3_MRI', '1_Anatomicals', '3D_Anat', strcat(Subject.ID, '_3Danat_brain.nii.gz'));
            flirt_output1  = fullfile(tmp_out_FilePath, append(Input2, '_to2Danat', '.nii.gz'));
            flirt_output2  = fullfile(tmp_out_FilePath, append(Input2, '_to3Danat', '.nii.gz'));
            flirt_init1    = fullfile(Paths.fMRIrunsPath, strcat(Subject.fMRIruns(1), '.feat'), 'reg', 'example_func2initial_highres.mat');
            flirt_init2    = fullfile(Paths.fMRIrunsPath, strcat(Subject.fMRIruns(1), '.feat'), 'reg', 'example_func2highres.mat');
            flirt_options  = '-paddingsize 0.0 -interp trilinear';
            
            % Run registration to 2Danat space
            command_flirt1 = ['flirt', '-in', flirt_input, '-applyxfm', '-init', flirt_init1, '-out', flirt_output1, flirt_options, '-ref', flirt_2Danat];
            command_flirt1 = strjoin(command_flirt1);
            system(command_flirt1);
            
            % Run registration to 3Danat space
            command_flirt2 = ['flirt', '-in', flirt_input, '-applyxfm', '-init', flirt_init2, '-out', flirt_output2, flirt_options, '-ref', flirt_3Danat];
            command_flirt2 = strjoin(command_flirt2);
            system(command_flirt2);
        end
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
