function [due_sig] = load_due_signal(path, sbj_code, task_str,fsamp)
% This script was created to assist students attending the course in:
% "Bioengineering - UNISS" in reading binary data collected with a commercially available device.
% 
%
% This script processes binary data collected
% with the commercially available device found in this link:
% https://www.otbioelettronica.it/index.php?option=com_content&view=article&id=96:due&catid=18:strumentazione-portatile&Itemid=101&lang=en
%
% Lines 27, 28, 29 may be changed according to which data users would like to process.
% This script specfically reads binary data and convert them to Voltage
% values as read from EMG and Auxiliary probes.  Data is stored in a Matlab
% structure named Probe.
%
% For example, data sampled by Channel 2 of Probe number 5 can be retrived
% with the commend "Probe(5).Ch2"
%
% For detailed information on the data collected with this device and on
% how to work with this data can be consulted in the Laboratory Guidelines.
% 
% Example on how to call this function
%
% Probe = ReadEMGAndAngleData('Data_201952914180',true,1,1); % in this file there is only electrogoniometer data
% Probe = ReadEMGAndAngleData('Data_20195291465',true,1,[1 2 5]); % in this file there are EMG and electrogoniometer data
%
% Author: Taian Vieira (tmartins@uniss.it)
% Last modified: 10/11/20 


%% Script starts here (do not change any lines from this point)

path_due= [path '\emg\'];
due_file= dir(fullfile(path_due, [sbj_code '*' task_str '_due*']));
filename= [path_due,due_file.name];
fid = fopen(filename,'r'); % open file for reading
data = fread(fid,[1 inf],'uint8'); % read binary data
fclose(fid); % close file

global CH_code

Param.ProbesConnected = 1:5;
Param.PcktSize = 265; % Bytes; dimension of single data packets
Param.PcktControl = 263; % Byte number; counter indicating whether packets have been lost
Param.ProbeByte = 262; % Byte number; byte within each packet with the number of the probe to which the packet belongs
Param.DataBytes = 5:260; % Byte numbers; data samples (see pdf material for how binary data are formated)
Param.EMGResolution = (2^16)-1; % Number of bits for the EMG probes (the most significant bit is reserved for signal
Param.Resolution = (2^16)-1; % Number of bits for the auxiliary probe (that with electrogoniometer data)
Param.DynRange = 3.3; % (V) AD Dynamic range
Param.Gain = 200; % (V/V) Gain for the EMG probes

%% Identification of data packets and extraction of physical data from each packet

indHead=find(data==255); % Identify headers

if isempty(indHead) % there is no data arrived
    return
end % if no packets were received

if numel(data)<indHead(end)+3 % check if all three bytes of header are present after the last occurrence of 255
    indHead(end) = []; % remove the last header
end

indHead255 = indHead(find(data(indHead+1)==255 & data(indHead+2)==0 & data(indHead+3)==0));% find indices starting sequences of 255 255 0 0

if numel(data)<(indHead255(end)+Param.PcktSize-1) % if there is not a full packet of data after the last Header
    data(indHead255(end):end) = []; % store the last, incomplete data packet to be read in the next iteration
    indHead255(end) = []; % remove the last header
end % process packets into the signals from single probes

if ~isempty(indHead255(diff(indHead255)~=265)) % if for any reason there is not full data packet between consecutive headers
    indHead255(diff(indHead255)~=265) = []; % remove the header leading to incomplete data packet
    samples = [[0:(Param.PcktSize-1)]' ones(Param.PcktSize,1)] * [ones(1,numel(indHead255)); indHead255]; % samples corresponding to single packets\
    data = data(samples);
else
    data = data(indHead255(1):(indHead255(end)+264));
    data = reshape(data,265,numel(indHead255));
end

%% looping across probes to extract data
for i_ch=1:length(CH_code)
    curr_ch=CH_code(i_ch);
    for i_muscle = Param.ProbesConnected
        
        %% getting data, converting it to mV unit, spliting into both channels and adding offset for visualization
        currentprobe = find(data(Param.ProbeByte,:)==i_muscle); % finding packets corresponding to the current probe
        
        %% next seven lines are necessary to identify any loss of packets for the current probe
        pcktscounter = data(Param.PcktControl,currentprobe); % counter used to identify missed packets
        pcktscounter = diff([pcktscounter(1)-1;pcktscounter(:)]); % pcktscounter should contain only 1 and -255, if there was no packet loss
        pcktscounter(pcktscounter<0)=pcktscounter(pcktscounter<0)+256; % replacing instance with negative values by 1, if no packet has been lost, or with the number of packets lost +1
        pcktscounter = cumsum(pcktscounter); % converting counter values to the actual number of processed packets
        Npckts = pcktscounter(end); % expected number of packets should none have been lost
        ProbeData = nan(numel(Param.DataBytes),Npckts); % initialising the matrix to store sampled data with "nans" (Not A Number)
        ProbeData(:,pcktscounter) = double(data(Param.DataBytes,currentprobe)); % replacing nan occurrences with packets actually read
        
        % ProbeData contains bytes with sampled data and NaN whenever there was a packet loss
        ProbeData = ProbeData(:); % reshaping data to a column vector
        ProbeData = (ProbeData(1:2:end,1)*256+ProbeData(2:2:end,1))*Param.DynRange/Param.Resolution; % converting data to V unit
        
        % assigning signals to each probe channel
        Ch = ProbeData(i_ch:2:end,1);
                
        %% next lines are necessary to filter EMGs and auxiliary data while considering for the occurrence of loss of packets
        nanind = find(isnan(Ch)); nonnanind = find(~isnan(Ch)); % indices pointing to NaN and data values respectivelly; NaN and data indices, are the same for both channels
        if ~isempty(nanind) % if there are instances of NaN values (i.e., if packets have been lost)
            Ch = interp1(nonnanind,Ch(nonnanind),1:numel(Ch),'pchip'); % interpolating NaN occurrences using picewise cubic method
            
        end
        Ch=(Ch-nanmean(Ch))/Param.Gain;
        
        Ch(nanind) = nan; % assigning NaN values to instances of packet loss
        due_sig(i_muscle).(curr_ch) = Ch; %non-filtrated signal
        
        
        N(i_muscle) = length(due_sig(i_muscle).(curr_ch)); % total number of samples that should have been collected for the current probe
    end
end
    
    % clear ProbeData currentprobe probes fid data indHead indHead255 samples pcktscounter Npckts
    
    %% The following commands are meant to plot data and to truncate data, ensuring the same number of samples to all channels
try
    N = min(N(N>0)); % smaller number of samples that should have been saved for a given probe
catch
    keyboard
end
    for i_ch = 1:length(CH_code)
        curr_ch=CH_code(i_ch);
        for i_muscle = 1:numel(Param.ProbesConnected)
            due_sig(Param.ProbesConnected(i_muscle)).(curr_ch) = due_sig(Param.ProbesConnected(i_muscle)).(curr_ch)(1:N); % truncating data to N samples
        end
    end
    
    %find the trigger transitions
    
    trig=find(due_sig(4).Ch2>0);
    due_sig(4)=[];
    
    %remove the signal epoch before and after the trigger transitions
    for i_ch=1:length(CH_code)
        curr_CH=CH_code(i_ch);
        for i_muscle=1:length(due_sig)
            due_sig(i_muscle).(curr_CH) = due_sig(i_muscle).(curr_CH)(trig(1):trig(end));
            
            if size(due_sig(i_muscle).(curr_CH),1)>size(due_sig(i_muscle).(curr_CH),2)
               due_sig(i_muscle).(curr_CH)=(due_sig(i_muscle).(curr_CH))';
            end
        end
    end
end

