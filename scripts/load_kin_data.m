%Load kinematic data
%
%Input:
%       fsamp: sampling frequency in samp/s
%Output:
%       tr: structure with one field fopr each marker
%       each field tr.(marker) contains a matrix [nsamps 3] 
%       with one line for each acquired frame 
%       Each line contains the x,y, and z coordinates of the marker 
%               
%               
function [tr]= load_kin_data(fsamp)

global markers

signal_subfolders= {'due_pro', 'meacs', 'vicon'};
task_list= {'gait','faster','bike_free','bike_stuck'};
selpath = uigetdir();
i= findstr(selpath,'\');
sbj_code= selpath(i(end)+1:end);

[indx,tf] = listdlg('PromptString',{'Select the task'},...
    'SelectionMode','single','ListString',task_list);
task_str= task_list{indx};

[traj]=load_Vicon_data(selpath, sbj_code, task_str);
traj= traj';
% global markers= {'LASI'; 'RASI'; 'LPSI'; 'RPSI';'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'; ...
%     'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'};
for i_m=1:length(markers)
    tr.(markers{i_m})= [traj((i_m-1)*3+1,:); traj((i_m-1)*3+2,:); traj((i_m-1)*3+3,:)]';
end

