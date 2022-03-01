% TESINA - Rigazio Sofia, Roccaro Lucia, Romano Anastasio, Ruzzante Elena
clear variables, close all, clc

addpath scripts

%% define global variables to reduce parameters in the functions
global markers
global segments
global joints
global sides
global Antro
global muscle_code % EMG data

% sampling frequencies
kin_fsamp = 100;
emg_fsamp = 2048;

%% Loading Data
markers = {'RASI'; 'LASI'; 'RPSI'; 'LPSI';'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'; ...
    'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'};
segments = {'pelvis'; 'femur'; 'thigh'; 'foot'};
joints = {'LH'; 'RH'; 'LK'; 'RK'; 'LA'; 'RA'};
coords = {'x'; 'y'; 'z'};
sides = {'left'; 'right'};

% Muscle codes
muscle_code = ["BF","GMed","TA","VL","SOLm", "GM"];
[sd_sig, traj] = load_kin_EMG_data();

% swap columns 1 and 2 in order to have medio lateral axis along y
for i_m = 1:length(markers)
    traj.(markers{i_m}) = traj.(markers{i_m})(:,[2 1 3]);
end

n_frames_kin = length(traj.LANK);
n_marker = length(markers);
n_samp_emg = length(sd_sig.BF);
n_muscles = length(muscle_code);

%% Definizione sistemi di riferimento
% Visualizzare la posizione dei marker in 3D per ogni frame acquisito
color = hsv(n_marker);  % per dare ad ogni marker un colore univoco

% dal lab02 sono stati quindi valutati i limiti degli assi
xmin = -4000; xmax = 6000; ymin = -2000; ymax = 200; zmin = 0; zmax = 1200;

%% Calcolo dei sistemi di riferimento delle articolazioni (orientamento e posizione)
% Antropometric masures (in mm) for hip joint center estimation;
Antro.LASI_RASI_dist = 240;  % distance between LASI and RASI
Antro.leg_length = 980;	    % leg length

% Antropometric masures (in mm) for knee joint center estimation
Antro.knee_width = 100;      % knee width in mm

% Antropometric masures (in mm) for ankle joint center estimation
Antro.ankle_width = 80;      % ankle width in mm

Antro.mDiameter = 10;  % marker diameter in mm

% attraverso la funzione 'calc_references.m' ricaviamo due struct con i
% sistemi di riferimento (orientamento e posizione)
[loc_ref, JC] = calc_references(traj);

%% STEP 1: Calcolare angoli articolari dell'anca, ginocchio e caviglia dell'arto sinistro
fprintf('STEP 1: calcolo angoli articolari di anca, ginocchio e caviglia dell''arto sinistro\n')
% i sistemi di riferimento calcolati dalla funzione calc_references.m e
% contenuti nella variabile loc_ref corrispondono alle coordinate della
% punta dei corrispondenti versori x,y,z rispetto all'origine (0,0,0).
% Si ottengono quindi 6 versori, tre del sistema prossimale e tre del
% distale, con origine coincidente.

% NOTA: abbiamo ruotato anche le direzioni del piede per fare in modo tale
% da avere tutte terne x,y,z destrorse con asse z verticale

angles = struct();
% 1: HIP - KNEE
% articolazione prossimale
coord_prox.x = loc_ref.hip.x;
coord_prox.y = loc_ref.hip.y;
coord_prox.z = loc_ref.hip.z;

% articolazione distale
coord_dist.x = loc_ref.knee_L.x;
coord_dist.y = loc_ref.knee_L.y;
coord_dist.z = loc_ref.knee_L.z;

% calcolo degli angoli articolari
[angles.hip.alpha, angles.hip.beta, angles.hip.gamma] = calc_angles(coord_prox, coord_dist);

% 2. KNEE - ANKLE
% articolazione prossimale
coord_prox.x = loc_ref.knee_L.x;
coord_prox.y = loc_ref.knee_L.y;
coord_prox.z = loc_ref.knee_L.z;

% articolazione distale
coord_dist.x = loc_ref.ankle_L.x;
coord_dist.y = loc_ref.ankle_L.y;
coord_dist.z = loc_ref.ankle_L.z;

% calcolo degli angoli articolari
[angles.knee_L.alpha, angles.knee_L.beta, angles.knee_L.gamma] = calc_angles(coord_prox, coord_dist);

% 3. ANKLE - FOOT (da chiedere)
% articolazione prossimale
coord_prox.x = loc_ref.ankle_L.x;
coord_prox.y = loc_ref.ankle_L.y;
coord_prox.z = loc_ref.ankle_L.z;

% articolazione distale
coord_dist.x = -loc_ref.foot_L.z;
coord_dist.y = loc_ref.foot_L.y;
coord_dist.z = loc_ref.foot_L.x;

% coord_dist.x = loc_ref.foot_L.x;
% coord_dist.y = loc_ref.foot_L.y;
% coord_dist.z = loc_ref.foot_L.z;

% calcolo degli angoli articolari
[angles.ankle_L.alpha, angles.ankle_L.beta, angles.ankle_L.gamma] = calc_angles(coord_prox, coord_dist);


%% STEP 2: Identificazione dei cicli del cammino utilizzando la traiettoria del marker sul tallone sinistro
fprintf('STEP 2: identificazione cicli del cammino\n')
% Ogni ciclo comincia e finisce con l’appoggio del tallone
asset = 0:1/kin_fsamp:(n_frames_kin-1)/kin_fsamp;
[~,locs] = findpeaks(-traj.LHEE(:,3), 'MinPeakProminence', 8); % la funzione findpeaks individua i massimi locali
duration = diff(asset(locs));

figure('Name','Istogramma delle durate di ogni ciclo', 'Position', [0, 300, 500, 400])
histogram(duration,30)
xline(mean(duration),'-.r', 'Media delle durate')
xlabel('time (s)'), ylabel('# cycles')
title('Istogramma delle durate di ogni ciclo')

n_cycles = 0;
i = 0;
while i < length(duration)
    i = i+1;
    if duration(i) < mean(duration) % escludiamo i cicli non completi (corrispondenti alle svolte)
        n_cycles = n_cycles +1;
        start_cycle(n_cycles) = locs(i); % samples
        end_cycle(n_cycles) = locs(i+1); % samples  
    else
        i = i+1; % scartiamo anche i primi cicli dopo la svolta
    end
end

% Visualizzazione in 3D in tutti i frame contemporaneamente
% figure('Name', '3D visualization')
% activation = traj.LHEE;
% plot3(activation(:,1), activation(:,2), activation(:,3)) % andamento del passo
% hold on
% scatter3(activation(start_cycle,1), activation(start_cycle,2),activation(start_cycle,3)) % istanti di inizio ciclo
% xlabel('x axis'), ylabel('y axis'), zlabel('z axis')
% hold on, axis equal
% grid on

% cerchiamo eventuali altri cicli non completi e li rimuoviamo
for i = 1:n_cycles
    trovato = 0; % questa variabile va a 1 se troviamo dei NaN -> ciclo non completo -> da eliminare
    for j = start_cycle(i):end_cycle(i)
        for k = 1:n_marker
            if (isnan(sum(traj.(markers{k})(j,:))) && trovato == 0)
                % se entriamo qui, abbiamo trovato dei NaN tra i marker
                start_cycle(i) = NaN;
                end_cycle(i) = NaN;
                trovato = 1;
            end
        end
    end
end

% rimuoviamo i NaN dai vettori start_cycles e end_cycles
start_cycle = start_cycle(~isnan(start_cycle));
end_cycle = end_cycle(~isnan(end_cycle));

% ricalcoliamo il numero di cicli
n_cycles = length(start_cycle);

% Visualizzazione della traiettoria del tallone -> marker LHEE con istanti
% di inzio e fine ciclo
figure('Name','Traiettoria LHEE','Position',[550, 300, 1000, 400])
plot(asset,traj.LHEE(:,3), 'b', 'DisplayName', 'LHEE')
hold on, grid on
p_end = plot(asset(end_cycle),traj.LHEE(end_cycle,3), '.r', 'MarkerSize', 20); % fine ciclo
p_start = plot(asset(start_cycle), traj.LHEE(start_cycle,3), '*b'); % inizio ciclo
xlabel('time (s)')
ylabel('z axis (mm)')
title('Traiettoria LHEE')
legend([p_start p_end],{'start','end'})

fprintf('Premere un tasto per continuare\n'), pause, clc
%% STEP 3: Calcolare l’andamento medio e la variabilità degli angoli articolari su tutti i passi
fprintf('STEP 3: calcolo andamento medio e variabilità angoli articolari\n')
% Creiamo una matrice contenente tutti i cicli del passo, ricampionati in
% modo che ogni ciclo abbia lo stesso numero di campioni

% Abbiamo osservato che il numero di campioni per ogni ciclo del passo
% varia tra 108 e 120. Abbiamo deciso quindi di sottocampionare tutti i
% cicli a 101 campioni, in modo da prendere in considerazione le
% percentuali del ciclo del passo da 0 a 100% in passi unitari.

cycle_length = 101; % campioni
% asse x per plottare in funzione della percentuale del ciclo del passo
asse_x = linspace(0,100,cycle_length);
% MEMO:
% - alpha = adduzione-abduzione
% - beta = flessione ed estensione
% - gamma = intra- ed extra-rotazione
for i = 1:n_cycles
    original_length = end_cycle(i)-start_cycle(i);

    % hip
    angles.hip.alpha_cycles(:,i) = resample(angles.hip.alpha(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);
    angles.hip.beta_cycles(:,i) = -resample(angles.hip.beta(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);
    angles.hip.gamma_cycles(:,i) = resample(angles.hip.gamma(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);

    % knee
    angles.knee_L.alpha_cycles(:,i) = resample(angles.knee_L.alpha(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);
    angles.knee_L.beta_cycles(:,i) = resample(angles.knee_L.beta(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);
    angles.knee_L.gamma_cycles(:,i) = resample(angles.knee_L.gamma(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);

    % ankle
    angles.ankle_L.alpha_cycles(:,i) = resample(angles.ankle_L.alpha(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);
    angles.ankle_L.beta_cycles(:,i) = -resample(angles.ankle_L.beta(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);
    angles.ankle_L.gamma_cycles(:,i) = resample(angles.ankle_L.gamma(start_cycle(i):end_cycle(i)-1),cycle_length,original_length);

end

% %% plot completo angoli
% figure
% for art = 1:3 % 3 articolazioni
%      switch art
%         case 1
%             joint = 'hip';
%         case 2
%             joint = 'knee_L';
%         case 3
%             joint = 'ankle_L';
%     end
%     for ang = 1:3 % 3 angoli
%         switch ang
%             case 1
%                 angolo = 'alpha';
%             case 2
%                 angolo = 'beta';
%             case 3
%                 angolo = 'gamma';
%         end
%         valued = eval(sprintf('angles.%s.%s_cycles',joint,angolo)) * 180/pi; % convertiamo in gradi
%         mean_valued = mean(valued,2);
%         st_err_valued = std(valued,[],2)/sqrt(size(valued,2));
% 
%         subplot(3,3,3*(art-1)+ang)
%         plot(asse_x,mean_valued,'k','LineWidth',2)
%         hold on, grid on
%         plot(asse_x,mean_valued + st_err_valued, 'r', LineStyle="--")
%         plot(asse_x,mean_valued - st_err_valued, 'r', LineStyle="--")
%         % barre di errore
%         x2 = [asse_x, fliplr(asse_x)];
%         inBetween = [(mean_valued'-st_err_valued'), fliplr(mean_valued'+st_err_valued')];
%         fill(x2, inBetween, 'r', 'FaceAlpha',0.3);
%         title(sprintf('%s - %s',joint, angolo))
%         hold off
%         xlabel('% del ciclo del passo')
%         ylabel('angolo (°)')
%     end
% end

fprintf('Premere un tasto per visualizzare gli angoli separatamente sui relativi piani\n'), pause, clc
%% Plot solo angoli sul piano sagittale
figure('Name', 'Flesso-estensione', 'Position', [10 500 1500 300])
sgtitle('Angoli sul piano sagittale')
for art = 1:3 % 3 articolazioni
        switch art
        case 1
            joint = 'hip';
            label = 'Hip Flexion/extension';
        case 2
            joint = 'knee_L';
            label = 'Knee Flexion/extension';
        case 3
            joint = 'ankle_L';
            label = 'Ankle Dorsi/Plantar flexion';
        end
    angolo = 'beta';
    valued = eval(sprintf('angles.%s.%s_cycles',joint,angolo)) * 180/pi; % convertiamo in gradi
    mean_valued = mean(valued,2);
    st_err_valued = std(valued,[],2)/sqrt(size(valued,2));

    subplot(1,3,art)
    plot(asse_x,mean_valued,'k','LineWidth',2)
    hold on, grid on
    plot(asse_x,mean_valued + st_err_valued, 'r', LineStyle="--")
    plot(asse_x,mean_valued - st_err_valued, 'r', LineStyle="--")
    % barre di errore
    x2 = [asse_x, fliplr(asse_x)];
    inBetween = [(mean_valued'-st_err_valued'), fliplr(mean_valued'+st_err_valued')];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.3);
    title(label)
    hold off
    xlabel('% del ciclo del passo')
    ylabel('angolo (°)')
end

%% Plot solo angoli piano trasverso
figure('Name', 'Piano trasverso', 'Position', [10 300 1500 300])
sgtitle('Angoli sul piano trasverso')
for art = 1:3 % 3 articolazioni
        switch art
        case 1
            joint = 'hip';
            label = 'Hip Intra-extra rotation';
        case 2
            joint = 'knee_L';
            label = 'Knee Intra-extra rotation';
        case 3
            joint = 'ankle_L';
            label = 'Ankle Intra-extra rotation';
        end
    angolo = 'gamma';
    valued = eval(sprintf('angles.%s.%s_cycles',joint,angolo)) * 180/pi; % convertiamo in gradi
    mean_valued = mean(valued,2);
    st_err_valued = std(valued,[],2)/sqrt(size(valued,2));

    subplot(1,3,art)
    plot(asse_x,mean_valued,'k','LineWidth',2)
    hold on, grid on
    plot(asse_x,mean_valued + st_err_valued, 'r', LineStyle="--")
    plot(asse_x,mean_valued - st_err_valued, 'r', LineStyle="--")
    % barre di errore
    x2 = [asse_x, fliplr(asse_x)];
    inBetween = [(mean_valued'-st_err_valued'), fliplr(mean_valued'+st_err_valued')];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.3);
    title(label)
    hold off
    xlabel('% del ciclo del passo')
    ylabel('angolo (°)')
end

%% Plot solo angoli piano frontale
figure('Name', 'intra-extra rotazione', 'Position', [10 100 1500 300])
sgtitle('Angoli sul piano frontale')
for art = 1:3 % 3 articolazioni
        switch art
        case 1
            joint = 'hip';
            label = 'Hip Adduction/Abduction';
        case 2
            joint = 'knee_L';
            label = 'Knee Varus/Valgus';
        case 3
            joint = 'ankle_L';
            label = 'Ankle Inversion/Eversion';
        end
    angolo = 'alpha';
    valued = eval(sprintf('angles.%s.%s_cycles',joint,angolo)) * 180/pi; % convertiamo in gradi
    mean_valued = mean(valued,2);
    st_err_valued = std(valued,[],2)/sqrt(size(valued,2));

    subplot(1,3,art)
    plot(asse_x,mean_valued,'k','LineWidth',2)
    hold on, grid on
    plot(asse_x,mean_valued + st_err_valued, 'r', LineStyle="--")
    plot(asse_x,mean_valued - st_err_valued, 'r', LineStyle="--")
    % barre di errore
    x2 = [asse_x, fliplr(asse_x)];
    inBetween = [(mean_valued'-st_err_valued'), fliplr(mean_valued'+st_err_valued')];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.3);
    title(label)
    hold off
    xlabel('% del ciclo del passo')
    ylabel('angolo (°)')
end
    
fprintf('Premere un tasto per continuare\n'), pause, clc

%% PLOT - Visualizzare la posizione e orientamento dei sistemi di riferimento
fprintf('Visualizzazione posizione e orientamento sistemi di riferimento\n')
% il vettore qui sotto ci permette invece nella visualizzazione di collegare correttamente i markers fra loro
link_markers = [1 2; 1 3; 2 4; 3 4; 5 6; 6 7; 7 8; 8 9; 9 10; 8 10; 11 12; 
    12 13; 13 14; 14 15; 15 16; 14 16];
scale = 100; % termine per scalare le frecce dei versori

% PLOT PRIMO FRAME CON SISTEMI DI RIFERIMENTO
figure('Position',[100, 100, 900, 800])
legend('AutoUpdate','off')
for j = 1:n_marker
    valued = traj.(markers{j});
    h_m(j) = plot3(traj.(markers{j})(500,1), traj.(markers{j})(500,2), traj.(markers{j})(500,3), '.k', 'MarkerSize', 10); 
    hold on, grid on
end

% bacino
plt_bacino = plot3(JC.hip(500,1), JC.hip(500,2), JC.hip(500,3), '.m', 'MarkerSize', 25, 'DisplayName', 'O_bacino');
plt_bacino_x = quiver3(JC.hip(500,1), JC.hip(500,2), JC.hip(500,3), scale*loc_ref.hip.x(500,1), scale*loc_ref.hip.x(500,2), scale*loc_ref.hip.x(500,3), 'r', 'DisplayName', 'asse x');
plt_bacino_y = quiver3(JC.hip(500,1), JC.hip(500,2), JC.hip(500,3), scale*loc_ref.hip.y(500,1), scale*loc_ref.hip.y(500,2), scale*loc_ref.hip.y(500,3), 'g', 'DisplayName', 'asse y');
plt_bacino_z = quiver3(JC.hip(500,1), JC.hip(500,2), JC.hip(500,3), scale*loc_ref.hip.z(500,1), scale*loc_ref.hip.z(500,2), scale*loc_ref.hip.z(500,3), 'b', 'DisplayName', 'asse z');

% ginocchio sinistro
plt_knee_L = plot3(JC.knee_L(500,1), JC.knee_L(500,2), JC.knee_L(500,3), '.c', 'MarkerSize', 25, 'DisplayName', 'O_ginocchio_sx');
plt_knee_L_x = quiver3(JC.knee_L(500,1), JC.knee_L(500,2), JC.knee_L(500,3), scale*loc_ref.knee_L.x(500,1), scale*loc_ref.knee_L.x(500,2), scale*loc_ref.knee_L.x(500,3), 'r');
plt_knee_L_y = quiver3(JC.knee_L(500,1), JC.knee_L(500,2), JC.knee_L(500,3), scale*loc_ref.knee_L.y(500,1), scale*loc_ref.knee_L.y(500,2), scale*loc_ref.knee_L.y(500,3), 'g');
plt_knee_L_z = quiver3(JC.knee_L(500,1), JC.knee_L(500,2), JC.knee_L(500,3), scale*loc_ref.knee_L.z(500,1), scale*loc_ref.knee_L.z(500,2), scale*loc_ref.knee_L.z(500,3), 'b');

% ginocchio destro
plt_knee_R = plot3(JC.knee_R(500,1), JC.knee_R(500,2), JC.knee_R(500,3), '.c', 'MarkerSize', 25, 'DisplayName', 'O_ginocchio_dx');
plt_knee_R_x = quiver3(JC.knee_R(500,1), JC.knee_R(500,2), JC.knee_R(500,3), scale*loc_ref.knee_R.x(500,1), scale*loc_ref.knee_R.x(500,2), scale*loc_ref.knee_R.x(500,3), 'r');
plt_knee_R_y = quiver3(JC.knee_R(500,1), JC.knee_R(500,2), JC.knee_R(500,3), scale*loc_ref.knee_R.y(500,1), scale*loc_ref.knee_R.y(500,2), scale*loc_ref.knee_R.y(500,3), 'g');
plt_knee_R_z = quiver3(JC.knee_R(500,1), JC.knee_R(500,2), JC.knee_R(500,3), scale*loc_ref.knee_R.z(500,1), scale*loc_ref.knee_R.z(500,2), scale*loc_ref.knee_R.z(500,3), 'b');

% caviglia sinistra
plt_ankle_L = plot3(JC.ankle_L(500,1), JC.ankle_L(500,2), JC.ankle_L(500,3), '.g', 'MarkerSize', 25, 'DisplayName', 'O_caviglia_sx');
plt_ankle_L_x = quiver3(JC.ankle_L(500,1), JC.ankle_L(500,2), JC.ankle_L(500,3), scale*loc_ref.ankle_L.x(500,1), scale*loc_ref.ankle_L.x(500,2), scale*loc_ref.ankle_L.x(500,3), 'r');
plt_ankle_L_y = quiver3(JC.ankle_L(500,1), JC.ankle_L(500,2), JC.ankle_L(500,3), scale*loc_ref.ankle_L.y(500,1), scale*loc_ref.ankle_L.y(500,2), scale*loc_ref.ankle_L.y(500,3), 'g');
plt_ankle_L_z = quiver3(JC.ankle_L(500,1), JC.ankle_L(500,2), JC.ankle_L(500,3), scale*loc_ref.ankle_L.z(500,1), scale*loc_ref.ankle_L.z(500,2), scale*loc_ref.ankle_L.z(500,3), 'b');

% caviglia destra
plt_ankle_R = plot3(JC.ankle_R(500,1), JC.ankle_R(500,2), JC.ankle_R(500,3), '.g', 'MarkerSize', 25, 'DisplayName', 'O_caviglia_dx');
plt_ankle_R_x = quiver3(JC.ankle_R(500,1), JC.ankle_R(500,2), JC.ankle_R(500,3), scale*loc_ref.ankle_R.x(500,1), scale*loc_ref.ankle_R.x(500,2), scale*loc_ref.ankle_R.x(500,3), 'r');
plt_ankle_R_y = quiver3(JC.ankle_R(500,1), JC.ankle_R(500,2), JC.ankle_R(500,3), scale*loc_ref.ankle_R.y(500,1), scale*loc_ref.ankle_R.y(500,2), scale*loc_ref.ankle_R.y(500,3), 'g');
plt_ankle_R_z = quiver3(JC.ankle_R(500,1), JC.ankle_R(500,2), JC.ankle_R(500,3), scale*loc_ref.ankle_R.z(500,1), scale*loc_ref.ankle_R.z(500,2), scale*loc_ref.ankle_R.z(500,3), 'b');

% piede sinistro
plt_toe_L = plot3(traj.LTOE(500,1), traj.LTOE(500,2), traj.LTOE(500,3), '.y', 'MarkerSize', 25, 'DisplayName', 'O_piede_sx');
plt_foot_L_x = quiver3(traj.LTOE(500,1), traj.LTOE(500,2), traj.LTOE(500,3), -scale*loc_ref.foot_L.z(500,1), -scale*loc_ref.foot_L.z(500,2), -scale*loc_ref.foot_L.z(500,3), 'r');
plt_foot_L_y = quiver3(traj.LTOE(500,1), traj.LTOE(500,2), traj.LTOE(500,3), scale*loc_ref.foot_L.y(500,1), scale*loc_ref.foot_L.y(500,2), scale*loc_ref.foot_L.y(500,3), 'g');
plt_foot_L_z = quiver3(traj.LTOE(500,1), traj.LTOE(500,2), traj.LTOE(500,3), scale*loc_ref.foot_L.x(500,1), scale*loc_ref.foot_L.x(500,2), scale*loc_ref.foot_L.x(500,3), 'b');

% piede destro
plt_toe_R = plot3(traj.RTOE(500,1), traj.RTOE(500,2), traj.RTOE(500,3), '.y', 'MarkerSize', 25, 'DisplayName', 'O_piede_dx');
plt_foot_R_x = quiver3(traj.RTOE(500,1), traj.RTOE(500,2), traj.RTOE(500,3), -scale*loc_ref.foot_R.z(500,1), -scale*loc_ref.foot_R.z(500,2), -scale*loc_ref.foot_R.z(500,3), 'r');
plt_foot_R_y = quiver3(traj.RTOE(500,1), traj.RTOE(500,2), traj.RTOE(500,3), scale*loc_ref.foot_R.y(500,1), scale*loc_ref.foot_R.y(500,2), scale*loc_ref.foot_R.y(500,3), 'g');
plt_foot_R_z = quiver3(traj.RTOE(500,1), traj.RTOE(500,2), traj.RTOE(500,3), scale*loc_ref.foot_R.x(500,1), scale*loc_ref.foot_R.x(500,2), scale*loc_ref.foot_R.x(500,3), 'b');

% % piede sinistro -> rotazione assi
% plt_foot_L_x = quiver3(traj.LTOE(500,1), traj.LTOE(500,2), traj.LTOE(500,3), -scale*loc_ref.foot_L.z(500,1), -scale*loc_ref.foot_L.z(500,2), -scale*loc_ref.foot_L.z(500,3), 'r');
% plt_foot_L_y = quiver3(traj.LTOE(500,1), traj.LTOE(500,2), traj.LTOE(500,3), scale*loc_ref.foot_L.y(500,1), scale*loc_ref.foot_L.y(500,2), scale*loc_ref.foot_L.y(500,3), 'g');
% plt_foot_L_z = quiver3(traj.LTOE(500,1), traj.LTOE(500,2), traj.LTOE(500,3), scale*loc_ref.foot_L.x(500,1), scale*loc_ref.foot_L.x(500,2), scale*loc_ref.foot_L.x(500,3), 'b');
% 
% % piede destro -> rotazione assi
% plt_foot_R_x = quiver3(traj.RTOE(500,1), traj.RTOE(500,2), traj.RTOE(500,3), -scale*loc_ref.foot_R.z(500,1), -scale*loc_ref.foot_R.z(500,2), -scale*loc_ref.foot_R.z(500,3), 'r');
% plt_foot_R_y = quiver3(traj.RTOE(500,1), traj.RTOE(500,2), traj.RTOE(500,3), scale*loc_ref.foot_R.y(500,1), scale*loc_ref.foot_R.y(500,2), scale*loc_ref.foot_R.y(500,3), 'g');
% plt_foot_R_z = quiver3(traj.RTOE(500,1), traj.RTOE(500,2), traj.RTOE(500,3), scale*loc_ref.foot_R.x(500,1), scale*loc_ref.foot_R.x(500,2), scale*loc_ref.foot_R.x(500,3), 'b');

legend([plt_bacino, plt_knee_L, plt_ankle_L, plt_toe_L, plt_bacino_x, plt_bacino_y, plt_bacino_z], ...
    {'origine bacino', 'origine ginocchio','origine caviglia', 'origine piede (TOE)','asse x', 'asse y', 'asse z'},'Location','southeast', 'AutoUpdate','off')
xlabel('x axis (mm)'), ylabel('y axis (mm)'), zlabel('z axis (mm)')
axis equal

for i_line = 1:size(link_markers,1)
        m_1 = link_markers(i_line,1);
        m_2 = link_markers(i_line,2);
    
        h_L(i_line) = line([traj.(markers{m_1})(500,1) traj.(markers{m_2})(500,1)], ...
            [traj.(markers{m_1})(500,2) traj.(markers{m_2})(500,2)], ...
            [traj.(markers{m_1})(500,3) traj.(markers{m_2})(500,3)],'Color','black','Linestyle','--');
end

% % scommentare per visualizzare sistema di riferimento plug-in-gait (vicon)
% % piede sinistro
% set(plt_foot_L_x, 'XData', traj.LTOE(500,1), 'YData', traj.LTOE(500,2), 'ZData', traj.LTOE(500,3), 'udata', scale*loc_ref.foot_L.x(500,1), 'vdata', scale*loc_ref.foot_L.x(500,2),'wdata', scale*loc_ref.foot_L.x(500,3))
% set(plt_foot_L_z, 'XData', traj.LTOE(500,1), 'YData', traj.LTOE(500,2), 'ZData', traj.LTOE(500,3), 'udata', scale*loc_ref.foot_L.z(500,1), 'vdata', scale*loc_ref.foot_L.z(500,2),'wdata', scale*loc_ref.foot_L.z(500,3))
% % piede destro
% set(plt_foot_R_x, 'XData', traj.RTOE(500,1), 'YData', traj.RTOE(500,2), 'ZData', traj.RTOE(500,3), 'udata', scale*loc_ref.foot_R.x(500,1), 'vdata', scale*loc_ref.foot_R.x(500,2),'wdata', scale*loc_ref.foot_R.x(500,3))
% set(plt_foot_R_z, 'XData', traj.RTOE(500,1), 'YData', traj.RTOE(500,2), 'ZData', traj.RTOE(500,3), 'udata', scale*loc_ref.foot_R.z(500,1), 'vdata', scale*loc_ref.foot_R.z(500,2),'wdata', scale*loc_ref.foot_R.z(500,3))

fprintf('Premere un tasto per continuare\n'), pause, clc
%% STEP 4: Calcolare l'andamento medio e la variabilità delle attivazioni muscolari
fprintf('STEP 4: calcolo andamento medio e variabilità attivazioni muscolari\n')
sig_f = calc_envelopes(sd_sig,emg_fsamp);

muscle_cycles = struct();

% riportiamo i vettori di inizio e fine dei cicli del passo alla frequenza
% di campionamento dei segnali emg
start_cycle_emg = round(start_cycle*emg_fsamp/kin_fsamp);
end_cycle_emg = round(end_cycle*emg_fsamp/kin_fsamp);

% Abbiamo osservato che il numero di campioni per ogni ciclo del passo
% varia tra 2355 e 2457. Abbiamo deciso quindi di sottocapionare tutti i
% cicli a 2300 campioni.

cycle_length = 2300; % campioni
for i = 1:n_cycles
    original_length = end_cycle_emg(i)-start_cycle_emg(i);

    % BF
    muscle_cycles.BF(:,i) = resample(sig_f.BF(start_cycle_emg(i):end_cycle_emg(i)-1),cycle_length,original_length)';
    % GMed
    muscle_cycles.GMed(:,i) = resample(sig_f.GMed(start_cycle_emg(i):end_cycle_emg(i)-1),cycle_length,original_length)';
    % TA
    muscle_cycles.TA(:,i) = resample(sig_f.TA(start_cycle_emg(i):end_cycle_emg(i)-1),cycle_length,original_length)';
    % VL
    muscle_cycles.VL(:,i) = resample(sig_f.VL(start_cycle_emg(i):end_cycle_emg(i)-1),cycle_length,original_length)';
    % SOLm
    muscle_cycles.SOLm(:,i) = resample(sig_f.SOLm(start_cycle_emg(i):end_cycle_emg(i)-1),cycle_length,original_length)';
    % GM_prox -> riga 1
    muscle_cycles.GMprox(:,i) = resample(sig_f.GM(1,start_cycle_emg(i):end_cycle_emg(i)-1),cycle_length,original_length)';
    % GM_dist -> riga 15
    muscle_cycles.GMdist(:,i) = resample(sig_f.GM(15,start_cycle_emg(i):end_cycle_emg(i)-1),cycle_length,original_length)';

end

% asse x per plottare in funzione della percentuale del ciclo del passo
asse_x = linspace(0,100,cycle_length);
signals = [muscle_code(1:5),'GMprox','GMdist'];

figure('Name', 'Attivazioni muscolari', 'Position', [100 100 1200 600])
sgtitle('Attivazioni muscolari')
for i = 1:7 % 7 segnali
    valued = eval(sprintf('muscle_cycles.%s',signals(i)));
    mean_valued = mean(valued,2);
    st_err_valued = std(valued,[],2)/sqrt(size(valued,2));

    subplot(4,2,i)
    plot(asse_x,mean_valued,'k','LineWidth',2)
    hold on
    plot(asse_x,mean_valued + st_err_valued, 'r', LineStyle="--")
    plot(asse_x,mean_valued - st_err_valued, 'r', LineStyle="--")
    % barre di errore
    x2 = [asse_x, fliplr(asse_x)];
    inBetween = [(mean_valued'-st_err_valued'), fliplr(mean_valued'+st_err_valued')];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.3);
    title(signals(i))
    hold off
    xlabel('% del ciclo del passo')
    ylabel('Ampiezza (V)')
end
