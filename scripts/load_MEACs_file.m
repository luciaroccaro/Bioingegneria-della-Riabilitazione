function signal= load_MEACs_file(handle,task_str,fsamp)

% %%
% %variables initialization: plot a gui (graphical user interface) in which is possible to set ADC
% %resolution [bit], ADC dynamic [V], Front-end gain [V/V] and Sampling
% %frequency [Hz]
% plotPrompt={'ADC resolution [bit]:','ADC dynamic [V]:','Front-end gain [V/V]:','Sampling frequency [Hz]:', 'Epoch length [ms]:', 'min Map Amp [uV_RMS]:', 'max Map Amp [uV_RMS]:', 'Interpolation? (0 - No, 1 - Yes)'};
% name='Plotting parameters';
% numlines=1;
% %Defaults values to fill the gui
% plotDefaultanswer={'16','2.4','192','2048', '500','0','50','1'};%set default answers inside text fields
% options.Resize='on';%make the just open window resizable
%
% %when the user press OK button, the values are recorded from the gui
% plotAnswer=inputdlg(plotPrompt,name,numlines,plotDefaultanswer,options);%put in answer the data just inserted from user

%parameters initialization from the gui
adcRes=16;
din=2.4;
gain=192;
fs=2048;
%epoch_length_ms = str2double(plotAnswer{5});
% min_rms = 0*1e-6;
% max_rms = 50*1e-6;
% interpFlag = 1;
%
% if interpFlag == 0
%     interpFlag = false;
% else
%     interpFlag = true;
% end

numChannels=33;
%calculate the number of ADC levels
maxLev=2^adcRes-1;
% epoch_length_samples = (epoch_length_ms/1000)*fs;


%to use for bluePlot wrote in C++ language
signal=fread(handle,[numChannels, inf], 'uint16',0,'l');   %to use with .sig files generated with bluePlotNem and the meacs system

for i=1:numChannels
    %converts the signal from level to Volts RTI (Referred To Input)
    signal(i,:)=((signal(i,:)/maxLev)*din)/gain;
    %eliminates the offset from the signal
    signal(i,:)=signal(i,:)-mean(signal(i,:));
    
end

% synchronized the signal using the trigger
trigger=signal(33,:);
if length(find(trigger>0))>length(find(trigger<0))
    s=find(trigger>0);
else
    s=find(trigger<0);
end
signal=signal(1:32,s(1):s(end));
end

