%% Section0: *** For this script to work, matlab needs to be opened through terminal
% e.g., in terminal, type /Applications/Matlab_R2019a.app/bin/matlab

% What it does - call on merge_func_runs function to merge together all functional runs
% to do this, example_func from Run2 onwards will be registered to the 
% example_func from run1. 
% In addition, the difference in timcourse means between that run and Run1
% subtracted from each voxel of the functional data for that run.
% once registration and normalization is complete, merging is done.

% Requirements
%
% Folder - Subject specific  (e.g., MinDch_Sub034)
%   Folder - one folder per run
%      filtered_func_ICA (i.e., one per run, in the run specific folder)
%      example_func (i.e., one per run, in the subject specific folder)
%   Folder - 'Merged_runs' that is currently empty
% 
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
tmp_Remove     = contains(SubjectFolders, '._');
tmp_keep       = find(tmp_Remove * -1 + 1);
SubjectFolders = SubjectFolders(tmp_keep);

% Generate GUI for the user to select the subjects they want to run
tmp_GUI_spacing = (30*numel(SubjectFolders)):-30:25;
tmp_GUI_height  = tmp_GUI_spacing(1) + 55;
tmp_GUI_spacing = tmp_GUI_spacing + 20;

tmp_fig               = uifigure;
tmp_fig.Visible       = 'on';
set(0,'units','pixels')
tmp_Pix_SS            = get(0,'screensize');
tmp_figCtr            = tmp_fig.Position(3:4)./2;
tmp_fig.Position(4)   = 600;
tmp_fig.Position(3)   = 350;
tmp_fig.Position(1:2) = (tmp_Pix_SS(3:4)./2) - tmp_figCtr;
tmp_fig.Name          = 'MIND Subject Selection';
tmp_fig.Scrollable    = 'on';
tmp_fig.Resize        = 'on';

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
    mkdir(fullfile(Paths.dataPath, SubjectFolders(CurrSub), '4_Analyses', 'Step1'))
    Subject.ID         = SubjectFolders(CurrSub);
    Paths.analysisPath = fullfile(Paths.dataPath, Subject.ID, '4_Analyses', 'Step1');

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
    

    % Make output folders
    Paths.registrationsPath = fullfile(Paths.analysisPath, 'Registrations');
    if ~isfolder(Paths.registrationsPath)
        mkdir(Paths.registrationsPath)
    end
    
    Paths.MergeFuncPath = fullfile(Paths.analysisPath, 'MergedFunctionals');
    if ~isfolder(Paths.MergeFuncPath)
        mkdir(Paths.MergeFuncPath)
    end
    
    clear tmp_*
    %% Section2: Merge Files
    
    % Merge the functional data together, then stop to allow time to check that this was done properly
    
    % modified by Tahereh, August 2024
    % Initialize and store all preprocessed functional fMRI files
    Subject.fMRIrunsICA = strings(1, numel(Subject.fMRIruns));

    for w = 1:numel(Subject.fMRIruns)
        line = fullfile(Paths.fMRIrunsPath, strcat(Subject.fMRIruns(w), '.feat'), strcat(Subject.ID, '_', Subject.fMRIruns(w), '_filtered_func_ICA.nii.gz'));
        Subject.fMRIrunsICA(w) = line;
    end

    clear 'w';

    
    Subject.fMRIrunsExampleFuncs = string(1:numel(Subject.fMRIruns))';
    for w = 1: numel(Subject.fMRIruns)
        line = fullfile(Paths.fMRIrunsPath, strcat(Subject.fMRIruns(w), '.feat'), 'reg', 'example_func.nii.gz');
        Subject.fMRIrunsExampleFuncs(w) = line;
    end
    
    Subject.fMRIrunsICAReg = string(1:numel(Subject.fMRIruns))';
    for r = 2:numel(Subject.fMRIruns)
        
        if isempty(Subject.fMRIruns(r)) == 0
            
        % This part computes the registration of between Run1 and all other
        % runs
        flirtToFirst_input = Subject.fMRIrunsExampleFuncs(r);
        flirtToFirst_ref = Subject.fMRIrunsExampleFuncs(1);
        flirtToFirst_out = fullfile(Paths.registrationsPath, strcat(Subject.ID, '_', Subject.fMRIruns(r), 'to', Subject.fMRIruns(1)));
        flirtToFirst_omat = fullfile(Paths.registrationsPath, strcat(Subject.ID, '_', Subject.fMRIruns(r), 'to', Subject.fMRIruns(1), '.mat'));
        flirtToFirst_options = '-bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear';
        
        command_flirttofirst = ['flirt', '-in', flirtToFirst_input, '-ref', flirtToFirst_ref, '-out', flirtToFirst_out, '-omat', flirtToFirst_omat, flirtToFirst_options];
        command_flirttofirst = strjoin(command_flirttofirst);
        system(command_flirttofirst);
        
        % This part applies the registrations
        regtofirst_in = strcat(Subject.fMRIrunsICA(r), ' -applyxfm');
        regtofirst_init = flirtToFirst_omat;
        regtofirst_out = fullfile(Paths.registrationsPath, strcat(Subject.ID, '_', Subject.fMRIruns(r), '_', 'filtered_func_ICA_to', Subject.fMRIruns(1), '.nii.gz'));
        regtofirst_options = '-paddingsize 0.0 -interp trilinear';
        regtofirst_ref = Subject.fMRIrunsExampleFuncs(1);
        
        command_regtofirst = ['flirt', '-in', regtofirst_in, '-init', regtofirst_init, '-out', regtofirst_out, regtofirst_options, '-ref', regtofirst_ref];
        command_regtofirst = strjoin(command_regtofirst);
        system(command_regtofirst);
        end
        
        Subject.fMRIrunsICAReg(r) = regtofirst_out;
    
    end
    
    % Normalize all runs to the first run
    %%% this part is going to calculate the timecourse mean (tmean) for all
    %%% runs, subtract the first Run tmean from the tmean of all other runs to find
    %%% a difference. 
    %%% The difference between the tmean of a run and Run1 will be subtracted
    %%% from that runs functional data (registered)
    
    %%% This calculates tmeans for Run1 (not registered) and all others
    %%% (registered)
    Subject.fMRIrunsTmean = string(1:numel(Subject.fMRIruns))';
    for run = 1:numel(Subject.fMRIruns)
        if run ==1 
        tmeans_in = Subject.fMRIrunsICA(run);
        tmeans_out = fullfile(Paths.MergeFuncPath, strcat(Subject.ID, '_', Subject.fMRIruns(run), '_tmean'));
        command_tmeans = ['fslmaths', tmeans_in, '-Tmean', tmeans_out];
        command_tmeans = strjoin(command_tmeans);
        system(command_tmeans);
        
        else
    
        tmeans_in = Subject.fMRIrunsICAReg(run);
        tmeans_out = fullfile(Paths.MergeFuncPath, strcat(Subject.ID, '_', Subject.fMRIruns(run), '_tmean'));
        command_tmeans = ['fslmaths', tmeans_in, '-Tmean', tmeans_out];
        command_tmeans = strjoin(command_tmeans);
        system(command_tmeans);
        
        end
        
        Subject.fMRIrunsTmean(run) = tmeans_out;
    
    end
    
    %%% This next part subtracts tmean for the first run from all other means.
    Subject.fMRIrunsTmeanDiff = string(1:numel(Subject.fMRIruns))';
    for run=2:numel(Subject.fMRIruns)
        
        tmeans_minusfirst_firstRun_in = Subject.fMRIrunsTmean(1);
        tmeans_minusfirst_nextRun_in = Subject.fMRIrunsTmean(run);
        tmeans_minusfirst_out = fullfile(Paths.MergeFuncPath, strcat(Subject.ID, '_', Subject.fMRIruns(run), '_tmean-', Subject.fMRIruns(1), '_tmean.nii.gz'));
        command_tmeans_minusfirst = ['fslmaths', tmeans_minusfirst_nextRun_in, '-sub', tmeans_minusfirst_firstRun_in, tmeans_minusfirst_out];
        command_tmeans_minusfirst = strjoin(command_tmeans_minusfirst);
        system(command_tmeans_minusfirst);
    
        Subject.fMRIrunsTmeanDiff(run) = tmeans_minusfirst_out;
    
    end
    
    %%% This part subtracts the difference in means from functional data for Runs2
    %%% and after
    Subject.fMRIrunsTmeanFinal = string(1:numel(Subject.fMRIruns))';
    for run=2:numel(Subject.fMRIruns)
    
        func_minus_tmeandif_in_func = Subject.fMRIrunsICAReg(run);
        func_minus_tmeandif_in_tmean = Subject.fMRIrunsTmeanDiff(run);
        func_minus_tmeandif_out = fullfile(Paths.MergeFuncPath, strcat(Subject.ID, '_', Subject.fMRIruns(run), '_filtered_func_ICA_reg-', Subject.fMRIruns(1), '_tmean', '.nii.gz'));
        command_func_minus_tmeandif = ['fslmaths', func_minus_tmeandif_in_func, '-sub', func_minus_tmeandif_in_tmean, func_minus_tmeandif_out];
        command_func_minus_tmeandif = strjoin(command_func_minus_tmeandif);
        system(command_func_minus_tmeandif);
    
        Subject.fMRIrunsTmeanFinal(run) = func_minus_tmeandif_out;
    
    end
    
    % Merge Run1 to all others
    %%% This takes 2 parts, one to merge run 1 and 2 and create a merged
    %%% file with merged filename. second to merge all subsequent runs to it. 
    
    %first, merge the first 2 runs, if there are that many.
    if numel(Subject.fMRIruns) > 1
        merge_first2_run1_in = Subject.fMRIrunsICA(1);
        merge_first2_run2_in = Subject.fMRIrunsTmeanFinal(2);
        merge_first2_out = fullfile(Paths.MergeFuncPath, strcat(Subject.ID, '_merged_functional_data_', Subject.fMRIruns(1), '-', Subject.fMRIruns(2),'.nii.gz'));
        command_merge12 = ['fslmerge -t', merge_first2_out, merge_first2_run1_in, merge_first2_run2_in];
        command_merge12 = strjoin(command_merge12);
        system(command_merge12);
    
        Subject.fMRImerged = merge_first2_out;
    % =====================
    % ADDED BY TAHEREH, APRIL 2024
    % we don't do the merging as there is only one fMRI run, we only create
    % the file using merge_first_2_out file to have the file with the name
    % that works for STEP2
    elseif numel(Subject.fMRIruns) == 1 % subject has only one fMIR run
        % Create an empty nifti file
        destin_dir = fullfile(Paths.dataPath, SubjectFolders(CurrSub), '3_MRI', '2_Functionals', '2_Preprocessed');
        empty_nifti_path = fullfile(destin_dir, 'empty.nii.gz');

        % Create a new header file with the desired pixel dimensions
        command_header_empty_file = ['fslcreatehd 64 64 24 1 3.750000 3.750000 5.000032 1 0 0 0 16 ' empty_nifti_path];
        system(char(command_header_empty_file));

        % Create an empty NIfTI file using the new header
        command_create_empty_file = ['fslmaths ' empty_nifti_path ' -mul 0 ' empty_nifti_path];
        system(char(command_create_empty_file));

        merge_first2_run1_in = Subject.fMRIrunsICA(1);
        merge_first2_run2_in = 1;
        merge_first2_out = fullfile(Paths.MergeFuncPath, strcat(Subject.ID, '_merged_functional_data_', Subject.fMRIruns(1), '-', Subject.fMRIruns(1),'.nii.gz'));     
        command_merge12 = ['fslmerge -t', merge_first2_out, merge_first2_run1_in, empty_nifti_path];
        command_merge12 = strjoin(command_merge12);
        system(command_merge12);
        % Remove the empty text file
        delete(empty_nifti_path);
        % =====================

    end
    
    % Second, merge the previously merged file with subsequent runs (if there
    % are more than two runs). Final filename ending in run1-run3 means that it
    % includes all runs from 1 to 3. 
    % this will keep looping through all functional runs, tacking on the
    % registered, normalized data to the end of the previously merged runs.
    
    if numel(Subject.fMRIruns) >2
        for run=3:numel(Subject.fMRIruns)
        merge_next_merged_in = Subject.fMRImerged;
        merge_next_next_in = Subject.fMRIrunsTmeanFinal(run);
        merge_next_out = fullfile(Paths.MergeFuncPath, strcat(Subject.ID, '_Merged_functional_data_', Subject.fMRIruns(1), '-', Subject.fMRIruns(run),'.nii.gz'));
        command_merge_next = ['fslmerge -t', merge_next_out, merge_next_merged_in, merge_next_next_in];
        command_merge_next = strjoin(command_merge_next);
        system(command_merge_next);
        
        end

    end

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