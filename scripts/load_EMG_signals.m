%LOad and preprocess  EMG signals recorded with both MEACS and DuePro systems.
%EMG acquired with DuePro are resampled to match the MEACS samples
function [sd_sig] = load_EMG_signals(path, sbj_code, task_str)

global CH_code

all_muscle_code=["BF","GMed","TA","VL","SOLm","SOLl","VMd","VMp", "GM"];

CH_code=["Ch1" "Ch2"];

%Due Signal
[due_sig]=load_due_signal(path, sbj_code, task_str);

n_probes= length(due_sig);
n_chs= length(CH_code);
for i_probe=1:n_probes
    for i_ch=1:n_chs
        curr_ch=CH_code(i_ch);
        sd_sig.(all_muscle_code((i_probe-1)*n_chs+i_ch))= due_sig(i_probe).(curr_ch);
    end
end
%Remove bad muscles
sd_sig= rmfield(sd_sig, "SOLl");
sd_sig= rmfield(sd_sig, "VMp");
sd_sig= rmfield(sd_sig, "VMd");


% Meacs Signal
mono_map= [1:2:32; 2:2:32]';

path_MEACs= [path '\emg\'];
MEACs_files= dir(fullfile(path_MEACs, [sbj_code '*' task_str '_meacs*.sig']));

for i_sig=1:length(MEACs_files)
    filePathName=[path_MEACs,MEACs_files(i_sig).name];
    fileHandle=fopen(filePathName,'r');

    if(findstr(filePathName, "M02BE4"))
        mono_MEACs.GM=load_MEACs_file(fileHandle,task_str);
    end
end

% total number of samples collected for the probes
N=[length(mono_MEACs.GM), length(due_sig(1).Ch1)];

N=min(N); % smaller number of samples that should have been saved for a given probe


mono_MEACs.GM =  mono_MEACs.GM(:,1:N); % truncating data to N samples
%Simulate 1x2 electrodes
for i_row=1:size(mono_map,1)
    signal(i_row,:)=mono_MEACs.GM(mono_map(i_row,1),:) + ...
        mono_MEACs.GM(mono_map(i_row,2),:);

    %Calc SD signals
    for i_row= 1:size(signal,1)-1
        sd_sig.GM(i_row,:)= signal(i_row+1,:)-signal(i_row,:);
    end
end