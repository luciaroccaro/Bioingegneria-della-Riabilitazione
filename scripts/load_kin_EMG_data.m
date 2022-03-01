%Read collected data
%Input: 
%   fsamp: sampling frequency in Hz
%Output: 
%   sd_sig
%   rhee_traj
function [sd_sig, tr]= load_kin_EMG_data()

global markers

signal_subfolders= {'emg', 'emg', 'vicon'};
task_list= {'gait'};
selpath = uigetdir();
i= findstr(selpath,'\');
sbj_code= selpath(i(end)+1:end);

[indx,tf] = listdlg('PromptString',{'Select the task'},...
    'SelectionMode','single','ListString',task_list);
task_str= task_list{indx};

[traj]=load_Vicon_data(selpath, sbj_code, task_str);
traj= traj';
for i_m=1:length(markers)
    tr.(markers{i_m})= [traj((i_m-1)*3+1,:); traj((i_m-1)*3+2,:); traj((i_m-1)*3+3,:)]';
end

%LOad EMG signals
[sd_sig] = load_EMG_signals(selpath, sbj_code, task_str);

end

