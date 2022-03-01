function sig_f = calc_envelopes(sd_sig,emg_fsamp)
plot_filtri = 0;

% struttura che conterrà i segnali processati
sig_f = struct();

%% interpolazione dei dati mancanti

sig_f.BF = fillmissing(sd_sig.BF,'linear');
sig_f.GMed = fillmissing(sd_sig.GMed,'linear');
sig_f.TA = fillmissing(sd_sig.TA,'linear');
sig_f.VL = fillmissing(sd_sig.VL,'linear');
sig_f.SOLm = fillmissing(sd_sig.SOLm,'linear');
sig_f.GM = fillmissing(sd_sig.GM,'linear',2);

%% filtraggio per ripulire i segnali grezzi

fNy = emg_fsamp/2; % Nyquist frequency

% filtro notch per eliminare la frequenza di rete -> filtro ricorsivo (ARMA)
z = 0.01;
B = 2;  % aumentando B si aumenta la profondità del rigetta banda
        % -> più B è alto più il segnale è pulito, ma essendo la banda
        %    maggiore si rischia di perdere parte del segnale
fc = 50;          % frequenza di rete
T = 1/emg_fsamp;  % periodo di campionamento
[b_notch, a_notch] = rico(z,B,fc,T);

% filtraggio (anticausale)
sig_f.BF = filtfilt(b_notch,a_notch,sig_f.BF);
sig_f.GMed = filtfilt(b_notch,a_notch,sig_f.GMed);
sig_f.TA = filtfilt(b_notch,a_notch,sig_f.TA);
sig_f.VL = filtfilt(b_notch,a_notch,sig_f.VL);
sig_f.SOLm = filtfilt(b_notch,a_notch,sig_f.SOLm);
sig_f.GM = filtfilt(b_notch,a_notch,sig_f.GM')';


% progettazione di un filtro passa-basso a 500Hz -> limite superiore banda
% del segnale EMG
Wp_low = 500/fNy;   % Passband frequency (normalized from 0 to 1)
Ws_low = 600/fNy;   % Stopband frequency (normalized from 0 to 1)
Rp_low = 1;         % Passband ripple (dB)
Rs_low = 20;        % Stopband attenuation (dB)
[Order_low,Wn_low] = buttord(Wp_low,Ws_low,Rp_low,Rs_low);
[b_low,a_low] = butter(Order_low,Wn_low);

% filtraggio (anticausale)
sig_f.BF = filtfilt(b_low,a_low,sig_f.BF);
sig_f.GMed = filtfilt(b_low,a_low,sig_f.GMed);
sig_f.TA = filtfilt(b_low,a_low,sig_f.TA);
sig_f.VL = filtfilt(b_low,a_low,sig_f.VL);
sig_f.SOLm = filtfilt(b_low,a_low,sig_f.SOLm);
sig_f.GM = filtfilt(b_low,a_low,sig_f.GM')';


% progettazione di un filtro passa-alto a 35Hz per eliminare i drift
Wp_high = 35/fNy;  % Passband frequency (normalized from 0 to 1)
Ws_high = 28/fNy;  % Stopband frequency (normalized from 0 to 1)
Rp_high = 1;       % Passband ripple (dB)
Rs_high = 5;       % Stopband attenuation (dB)
[Order_high,Wn_high] = buttord(Wp_high,Ws_high,Rp_high,Rs_high);
[b_high,a_high] = butter(Order_high,Wn_high, 'high');

% filtraggio (anticausale)
sig_f.BF = filtfilt(b_high,a_high,sig_f.BF);
sig_f.GMed = filtfilt(b_high,a_high,sig_f.GMed);
sig_f.TA = filtfilt(b_high,a_high,sig_f.TA);
sig_f.VL = filtfilt(b_high,a_high,sig_f.VL);
sig_f.SOLm = filtfilt(b_high,a_high,sig_f.SOLm);
sig_f.GM = filtfilt(b_high,a_high,sig_f.GM')';


if plot_filtri == 1
    figure
    freqz(b_notch,a_notch,emg_fsamp,emg_fsamp)
    title('filtro notch - 50Hz')

    figure
    freqz(b_low,a_low,emg_fsamp,emg_fsamp)
    title(sprintf('Low pass filter - cutoff frequency = %d Hz',Wp_low*fNy))

    figure
    freqz(b_high,a_high,emg_fsamp,emg_fsamp)
    title(sprintf('High pass filter - cutoff frequency = %d Hz',Wp_high*fNy))

end

%% estrazione dell'inviluppo

% rettifichiamo il segnale
sig_f.BF = abs(sig_f.BF);
sig_f.GMed = abs(sig_f.GMed);
sig_f.TA = abs(sig_f.TA);
sig_f.VL = abs(sig_f.VL);
sig_f.SOLm = abs(sig_f.SOLm);
sig_f.GM = abs(sig_f.GM);


% progettazione di un filtro passa-basso a 5 Hz
Wp_low = 5/fNy;     % Passband frequency (normalized from 0 to 1)
Ws_low = 9/fNy;     % Stopband frequency (normalized from 0 to 1)
Rp_low = 1;         % Passband ripple (dB)
Rs_low = 15;        % Stopband attenuation (dB)
[Order_low_1,Wn_low_1] = buttord(Wp_low,Ws_low,Rp_low,Rs_low);
[b_low_1,a_low_1] = butter(Order_low_1,Wn_low_1);

% progettazione di un filtro passa-basso a 3.8 Hz
Wp_low = 3.8/fNy;   % Passband frequency (normalized from 0 to 1)
Ws_low = 6/fNy;     % Stopband frequency (normalized from 0 to 1)
Rp_low = 1;         % Passband ripple (dB)
Rs_low = 15;        % Stopband attenuation (dB)
[Order_low_2,Wn_low_2] = buttord(Wp_low,Ws_low,Rp_low,Rs_low);
[b_low_2,a_low_2] = butter(Order_low_2,Wn_low_2);

% filtraggio (anticausale) con filtro LPF a 3.8 Hz
sig_f.BF = filtfilt(b_low_2,a_low_2,sig_f.BF);
sig_f.TA = filtfilt(b_low_2,a_low_2,sig_f.TA); 

% filtraggio (anticausale) con filtro LPF a 5Hz
sig_f.GMed = filtfilt(b_low_1,a_low_1,sig_f.GMed);
sig_f.VL = filtfilt(b_low_1,a_low_1,sig_f.VL);  
sig_f.SOLm = filtfilt(b_low_1,a_low_1,sig_f.SOLm);
sig_f.GM = filtfilt(b_low_1,a_low_1,sig_f.GM')';

if plot_filtri == 1
    figure
    freqz(b_low_1,a_low_1,emg_fsamp,emg_fsamp)
    title(sprintf('Low pass filter - cutoff frequency = %d Hz',Wp_low*fNy))

    figure
    freqz(b_low_2,a_low_2,emg_fsamp,emg_fsamp)
    title(sprintf('Low pass filter - cutoff frequency = %.1f Hz',Wp_low*fNy))
end

end