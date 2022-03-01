function [traj]=load_Vicon_data(selpath, sbj_code, task_str)

% %select the trajectory of the heel marker and resample using due signal
path_vicon = [selpath '\vicon_data\'];
filename= [sbj_code '_' task_str '_traj']; %load markers trajectories
load(fullfile(path_vicon,filename))

% filename= [sbj_code '_' task_str '_joint']; %load joint angle
% load(fullfile(path_vicon,filename))

%selected knee flec/ext, hip flec/ext and ankle fle/ext and resample with
%due signal frequency = 1000


end