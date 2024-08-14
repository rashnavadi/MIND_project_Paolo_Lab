# MIND_project_Paolo_Lab
modified versions of the five steps scripts for EEG-fMRI analysis on epilepsy patients

This MATLAB script is designed to automate the processing of fMRI and EEG data for subjects in a project, "MIND" study. Here’s a summary of what the script does:

Step 1 script:
1. Initial Setup
•	Paths and Directories:
o	Establishes paths based on the location of the script.
o	Warns the user to run MATLAB from the terminal to ensure the script functions correctly.
•	Subject Selection:
o	Searches for all potential subject folders within a specified directory.
o	Filters these folders to include only those relevant to the "MIND" study.
•	GUI for Subject Selection:
o	Generates a graphical user interface (GUI) allowing the user to select which subjects to include in the analysis.
o	The user selects subjects, and the script continues based on their choices.
2. Processing Each Subject
•	For each selected subject:
o	Directory Creation:
	Creates directories needed for storing analysis outputs.
o	fMRI and EEG Data Handling:
	Identifies and filters fMRI and EEG runs for each subject, ignoring certain files based on naming conventions.
o	Event Data Handling:
	Identifies event files relevant to the subject's EEG data.
o	Output Folder Creation:
	Creates additional directories needed for registrations and merged functional data.
3. Merging fMRI Data
•	fMRI Run Processing:
o	Handles cases where the subject has only one fMRI run (special case processing).
o	Registers and normalizes fMRI data from different runs relative to the first run.
•	Mean Calculation and Normalization:
o	Calculates the mean timecourse (tmean) for each fMRI run.
o	Normalizes all fMRI runs to the first run by subtracting the mean difference.
•	Merging Runs:
o	Merges all fMRI runs into a single file, normalizing and registering them to the first run.
o	If only one run is present, it creates an appropriate file for later steps.
4. Completion
•	Displays a message indicating the completion of processing for each subject.
5. Minor Functions
•	AddSubject Function:
o	Used within the GUI to update the list of selected subjects based on user input.
Step 2 script:
This MATLAB script is a follow-up to the previous one and is primarily designed to process EEG event data in relation to fMRI data. Here's a breakdown of what this script does:
Section 0: Introduction
•	Objective:
o	The script aims to merge EEG event files across different fMRI runs, adjusting the timing of events so that they are in "absolute time." This means that the events from subsequent runs are adjusted so that they don't start at time 0, but at the time when the previous run ended.
o	Multiple text files will be generated, each containing an increasing number of events.
Section 1: Initialization and GUI for Subject Selection
•	Paths Setup:
o	The script starts by setting up the paths needed for the project, similar to the previous script.
o	A warning message is displayed to ensure that MATLAB is run from the terminal.
•	Subject Selection:
o	The script identifies all relevant subject folders and generates a GUI for the user to select which subjects to process. This part is very similar to the first script.
Processing Each Subject
•	For each selected subject:
o	Path and Directory Setup:
	Directories specific to "Step 2" of the analysis are created.
	Paths for fMRI runs, EEG runs, and event files are set up.
o	Event File Processing:
	For each event type associated with a subject:
	Data Loading:
	EEG data for each run is loaded, along with the corresponding event timing files.
	The script processes the EEG data, including calculating epochs around the events and computing mean and standard deviation of the signal.
	Bipolar Channels:
	The script generates bipolar EEG channels (subtracting adjacent channels) and relabels them.
	A GUI is used to allow the user to select which bipolar channel to evaluate.
	Epoching and Statistics:
	The script epochs the EEG data around each event, then calculates and normalizes the mean and standard deviation of the signal.
	It also generates plots of the average waveform, including the mean and standard deviation, and saves these plots as PNG files.
Section 3: Event Timing Adjustment and Merging
•	Merging Event Timings:
o	The script reads event timings from the text files for each run and adjusts the timings so that they are continuous across runs (i.e., the time in run 2 does not reset to 0 but continues from the end of run 1).
o	The merged event timings are then saved into a new text file.
•	Volume Adjustment:
o	The script checks the total number of volumes in the fMRI data to ensure that the event timings do not extend beyond the length of the fMRI data.
o	Events occurring after the end of the fMRI data are removed.
•	Event Binning:
o	The script calculates how many events should be included in each output file, depending on the total number of events.
o	These event timings are then split into different files, each containing an increasing number of events, and saved as separate text files.
o	The timing of the last event in each file is also logged.
•	Final Output:
o	The merged event timings and the timing logs are saved to the appropriate directories.
•	Clearing Variables:
o	After processing each event type, variables are cleared to ensure that there is no carryover between subjects.
Minor Functions
•	AddSubject Function:
o	This function is used within the GUI to update the list of selected subjects based on user input.
Key Points:
•	Event Time Adjustment: The script ensures that events are aligned in absolute time across different runs.
•	Data Integrity Checks: The script includes various checks to handle cases where there might be only one run, no events, or incomplete data.
•	Output Generation: It creates multiple text files, each with an increasing number of events, and logs the timings of the last event in each file.
This script is designed to automate the process of merging and adjusting EEG event timings relative to fMRI data, ensuring that the data is correctly aligned for further analysis.
Step 3 script:
Overview of Step 3 Script
This script is designed to trim the merged fMRI data according to the lengths required for different increments of EEG events, based on the timings generated by previous scripts. The script uses FSL tools to perform the trimming and prepares the data for further analysis.
Detailed Breakdown:
Section 0: Introduction
•	Objective:
o	The script calls on FSL’s fslroi to trim the merged fMRI data according to the durations needed for various increments of EEG events.
o	It only runs if MATLAB is launched from the terminal, ensuring that the environment is set up correctly for FSL commands.
•	Requirements:
o	Subject-Specific Folder Structure:
	Merged_runs Folder: Contains the merged fMRI data and a .mat file with the run names.
	Event Timing Folder: Contains files with the timing of the last event and the total number of events for different increments.
Section 1: Initialization and GUI for Subject Selection
•	Paths Setup:
o	The script initializes paths based on the script’s location and sets up the data paths needed for the analysis.
o	A warning message is displayed, ensuring that the user runs MATLAB from the terminal.
•	Subject Selection:
o	The script identifies all relevant subject folders and generates a GUI for the user to select which subjects to process.
Processing Each Subject
•	Path and Directory Setup:
o	The script creates directories for "Step 3" of the analysis and sets up paths for fMRI runs, EEG runs, and event files.
•	Handling Single or Multiple fMRI Runs:
o	The script checks if the merged fMRI data file exists. If it doesn’t, it assumes there’s only one fMRI run and uses the corresponding preprocessed file instead.
Section 2: Trimming Functional Runs
•	Loading Event Timing Data:
o	The script loads the last event timings and the number of events per iteration, which were generated in the previous steps.
•	Trimming fMRI Data:
o	For each event increment:
	The script calculates the number of volumes needed based on the timing of the last event in that increment. This is done by dividing the last event timing by the TR (1.5 seconds) and rounding up.
	It then uses the fslroi command to trim the fMRI data, keeping only the volumes needed for that specific event increment. The trimmed data is saved as a new NIfTI file.
•	Saving Trimmed Volumes Information:
o	The script saves the number of volumes used for each iteration in a text file, which will be used in the next analysis step.
Summary of Key Steps:
1.	Initialization and GUI for Subject Selection: Sets up paths and allows the user to select subjects for processing.
2.	Handling Single or Multiple fMRI Runs: Ensures that the correct fMRI data is used, whether there’s one or multiple runs.
3.	Trimming fMRI Data: Uses FSL to trim the fMRI data according to the event timings, creating new files for each event increment.
4.	Saving Outputs: Saves the number of volumes used for each trimmed data file, preparing the data for subsequent analysis steps.
Purpose of the Script:
The purpose of this step is to ensure that the fMRI data is trimmed to precisely match the timing of specific EEG events. This allows for more accurate analysis in later steps, as the fMRI data will be aligned with the specific periods of interest defined by the EEG event timings. The script automates the process of generating these trimmed fMRI datasets, making the analysis process more efficient and reliable.

Step 4 script:
The Step 4 script is designed to create and execute FSL design files (.fsf) for fMRI analysis, followed by the amalgamation of different hemodynamic response functions (HRFs). Here’s a detailed explanation of what the script does and where you might need to adjust it to ensure full brain coverage in the output images:
Overview of Script Functions
1.	Initial Setup (Section0 & Section1):
o	The script initializes paths and checks if MATLAB is being run from the terminal.
o	It gathers a list of subject folders and prepares for user selection of subjects to process.
o	It sets up directories and paths for storing analysis outputs.
2.	Design File Creation (Section2):
o	For each subject and event, the script creates .fsf files. These files define how FSL’s feat will process the fMRI data.
o	It loads in templates and modifies them to match the specific number of volumes and events for each iteration.
o	The script prepares the output paths for storing results.
3.	Run fMRI Analysis (Section3):
o	The script loops through the generated .fsf files and runs FSL's feat to execute the analysis on the fMRI data.
o	This section essentially processes the fMRI data according to the design specifications created in Section2.
4.	Amalgamation and Registration (Section4):
o	After running the fMRI analysis, the script amalgamates different HRF models (referred to as "flobs").
o	It reads in the thresholded z-statistics images (e.g., thresh_zstat1.nii.gz) and creates composite images by taking the maximum value across the different HRF models.
o	The composite images are registered to both 2D and 3D anatomical spaces using FSL's flirt.
 
Step 5 script:
The Step 5 script is a final stage in analysis pipeline, focusing on comparing the results from different event iterations with the full event analysis using statistical methods. Here’s what the script does, broken down by its sections:
Overview of Script Functions
1.	Initial Setup (Section0 & Section1):
o	The script initializes paths and verifies that MATLAB is being run from the terminal.
o	It lists the subject folders available for analysis and allows the user to select which subjects to process.
o	Paths are set up for storing the outputs, including the analysis results, screenshots, and statistical tables.
2.	Running Statistics and Generating Outputs (Section2):
o	The script compares the fMRI analysis results from different event iterations with the analysis that includes all events.
o	This comparison is based on clusters identified in the fMRI data (e.g., significant activations or regions of interest).
Specific Steps in this Section:
o	Cluster Binarization:
	The clusters from the full event analysis are binarized. In this context, binarization means that voxels above a certain threshold are set to a value (e.g., 2), indicating they are part of a cluster.
	This binarized map is used as a reference to compare against the results from smaller iterations.
o	Loop through Event Iterations:
	For each event iteration (except the full event iteration), the script:
	Generates a binarized cluster map similar to the full event map but based on the current iteration.
	Takes screenshots of the cluster maps overlaid on the anatomical image using fsleyes. These screenshots help visualize how the brain activation patterns compare across different numbers of events.
	Calculates a confusion matrix and several statistical measures (e.g., True Positives, False Positives, Sensitivity, Specificity, F1 Score). These metrics quantify how well the iteration results match the full event analysis.
o	Confusion Matrix and Statistical Measures:
	The confusion matrix categorizes voxels as true positives, false positives, true negatives, or false negatives based on the comparison between the iteration and full event maps.
	Various metrics (e.g., Accuracy, Sensitivity, Specificity, F1 Score) are calculated to assess the performance of each iteration compared to the full event analysis.
o	Saving the Results:
	Screenshots for each iteration are saved to visualize the clusters.
	A statistical table is created, compiling all the metrics for each iteration. This table is saved as an Excel file, allowing further analysis or reporting.
3.	Minor Functions:
o	The script includes a minor function to manage the user interface for selecting subjects.
4.	One of the main changes in this step:
%         screenshot_command = append("fsleyes render --outfile ", outfile, " -s lightbox ", anat_filename, " ",  iteration_cluster_map_path_Screenshot, " -zx 3 --sliceSpacing 8 --zrange 5 100 --ncols 5 --nrows 3 -cm red-yellow");
% for fsleyes versioin like FSLeyes version 1.12.4 and probably higher versions the line below works
screenshot_command = append("fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 0.1 --zrange 0 1 --ncols 5 --nrows 3 --outfile ", outfile, " ", anat_filename, " ", iteration_cluster_map_path_Screenshot, " -cm red-yellow");
Key Features and Output:
•	Comparative Analysis: The main purpose of this step is to compare the analysis results from various iterations (different numbers of events) with the analysis including all events. It helps in understanding how many events are needed to achieve results comparable to using all events.
•	Statistical Evaluation: The confusion matrix and subsequent statistics provide a detailed assessment of how well each iteration captures the significant brain activity as identified by the full event analysis.
•	Visualization: The screenshots generated show the spatial distribution of brain activity for each iteration, overlaid on anatomical images. This visual output is crucial for quickly assessing whether the iterations are capturing the full brain activity or missing important areas.
•	Final Report: The Excel file generated at the end consolidates all the statistical metrics, making it easy to identify which event iterations are most accurate compared to the full event analysis.
           

![image](https://github.com/user-attachments/assets/e58e12d5-30ef-415d-a8c5-21dbce968b6e)
