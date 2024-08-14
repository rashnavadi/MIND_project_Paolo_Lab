MIND Study - fMRI and EEG Data Processing Scripts
This repository contains MATLAB scripts designed to automate the processing of fMRI and EEG data for subjects in the "MIND" study. The scripts guide users through multiple steps, from initial setup to final analysis, ensuring that data is prepared and processed correctly for further study.

Step 1: Initial Setup and Data Processing
Purpose: Automates the initial setup and processing of fMRI and EEG data for selected subjects.

Key Features:
Path Setup: Automatically establishes necessary paths and directories.
Subject Selection: GUI-based selection of subjects relevant to the "MIND" study.
Data Handling: Identifies, filters, and processes fMRI and EEG data, including event file management.
Output Creation: Generates directories and files needed for further analysis.
Processing Details:
Registers and normalizes fMRI data across runs.
Merges fMRI runs, adjusting them relative to the first run.
Outputs processed data for each subject with appropriate notifications.

Step 2: EEG Event Data Processing
Purpose: Processes EEG event data in relation to fMRI runs, ensuring correct timing alignment.

Key Features:
Event Timing Adjustment: Merges EEG event files across fMRI runs, adjusting timings to absolute time.
Event File Processing: Loads, processes, and analyzes EEG data, generating plots and statistics.
Output Generation: Creates multiple text files with increasing numbers of events, ensuring proper alignment with fMRI data.
Processing Details:
Adjusts and merges event timings across runs.
Ensures data integrity by handling cases with incomplete or single-run data.
Produces output files with adjusted event timings for further analysis.

Step 3: fMRI Data Trimming
Purpose: Trims merged fMRI data according to EEG event timings for precise analysis.

Key Features:
Data Trimming: Uses FSL tools to trim fMRI data according to the durations needed for EEG events.
Single/Multi-Run Handling: Checks for merged fMRI data and handles single/multiple runs accordingly.
Output Files: Generates new NIfTI files for each event increment, ready for subsequent analysis.
Processing Details:
Calculates the number of volumes needed based on event timings.
Trims and saves fMRI data, preparing it for specific event-related analysis.

Step 4: fMRI Analysis Design and Amalgamation
Purpose: Creates and executes FSL design files for fMRI analysis, followed by HRF amalgamation.

Key Features:
Design File Creation: Generates .fsf files for each subject and event, tailored to specific volumes and events.
fMRI Analysis Execution: Runs FSL's feat for fMRI data processing based on design specifications.
HRF Amalgamation: Combines different hemodynamic response functions, registering composite images to anatomical spaces.
Processing Details:
Automates the creation and execution of fMRI analysis designs.
Registers and amalgamates HRFs for comprehensive analysis.

Step 5: Comparative Analysis and Statistical Evaluation
Purpose: Compares results from different event iterations with full event analysis using statistical methods.

Key Features:
Cluster Analysis: Binarizes clusters from the full event analysis for comparison.
Event Iteration Comparison: Evaluates how well each iteration matches the full event analysis.
Statistical Outputs: Generates confusion matrices, accuracy, sensitivity, specificity, F1 scores, and visual comparisons.
Processing Details:
Visualizes brain activity distributions and compares them across different event counts.
Provides a detailed statistical report, summarizing the performance of each iteration.
