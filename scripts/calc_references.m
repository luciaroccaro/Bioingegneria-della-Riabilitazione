function [loc_ref, JC] = calc_references(traj) 
    global Antro
    JC = struct();
    loc_ref = struct();

    %% HIP
    LASI = traj.LASI;
    RASI = traj.RASI;
    RPSI = traj.RPSI;
    LPSI = traj.LPSI;
    
    % calcolo coordinate punto medio tra LASI e RASI
    mASI = (LASI+RASI)./2;
    
    % versore da RASI a LASI -> versore medio-laterale -> "asse y"
    v_med_lat = (LASI-RASI)./sqrt(sum(((LASI-RASI).^2),2));
    
    % calcolo coordinate punto medio tra RPSI e LPSI
    mPSI = (RPSI+LPSI)./2;
    
    % calcolo versore d'appoggio che va da mASI a mPSI
    v_mASI_mPSI = (mPSI-mASI)./sqrt(sum(((mPSI-mASI).^2),2));
    
    % prodotto vettoriale -> versore infero-superiore -> "asse z"
    v_inf_sup = cross(v_med_lat, v_mASI_mPSI);
    
    % versore postero-anteriore -> "asse x"
    v_post_ant = cross(v_med_lat, v_inf_sup);
    
    % salvataggio delle coordinate
    loc_ref.hip.x = v_post_ant;
    loc_ref.hip.y = v_med_lat;
    loc_ref.hip.z = v_inf_sup;
    
    theta = 0.5;    % in radianti
    beta = 0.314;   % in radianti
    Antro.Asis_troc_distance = 0.1288*Antro.leg_length - 48.56;
    MeanLegLength = Antro.leg_length; % il valore medio della lunghezza delle gambe lo consideriamo come la lunghezza della gamba stessa
    C = MeanLegLength*0.115-15.3;
    aa = Antro.LASI_RASI_dist/2;
    
    LHJC(1) = C*cos(theta)*sin(beta)-(Antro.Asis_troc_distance+Antro.mDiameter)*cos(beta);
    LHJC(2) = -(C*sin(theta)-aa);
    LHJC(3) = -C*cos(theta)*cos(beta)-(Antro.Asis_troc_distance+Antro.mDiameter)*sin(beta);
    LHJC = LHJC + mASI;

    RHJC(1) = C*cos(theta)*sin(beta)-(Antro.Asis_troc_distance+Antro.mDiameter)*cos(beta);
    RHJC(2) = (C*sin(theta)-aa);
    RHJC(3) = -C*cos(theta)*cos(beta)-(Antro.Asis_troc_distance+Antro.mDiameter)*sin(beta);
    RHJC = RHJC + mASI;

    % origine bacino -> punto medio delle coordinate RHJC e LHJC
    JC.hip = [(LHJC(:,1)+RHJC(:,1))/2, (LHJC(:,2)+RHJC(:,2))/2, (LHJC(:,3)+RHJC(:,3))/2];

    %% KNEE
    KneeOS = (Antro.mDiameter+Antro.knee_width)/2;
    
    % +++++ Left knee +++++++++++++++++++++++++++++++++
    v_L_THI_KNE = (traj.LTHI-traj.LKNE)./sqrt(sum(((traj.LTHI-traj.LKNE).^2),2));
    v_L_HJC_KNE = (LHJC-traj.LKNE)./sqrt(sum(((LHJC-traj.LKNE).^2),2));
    
    v_post_ant = cross(v_L_THI_KNE, v_L_HJC_KNE)./sqrt(sum(((v_L_HJC_KNE-v_L_THI_KNE).^2),2));

    JC.knee_L = chordPiG(traj.LTHI,LHJC,traj.LKNE,KneeOS);
    
    %versore che va da LHJC a LKJC (asse z)
    v_inf_sup = (LHJC-JC.knee_L)./sqrt(sum(((LHJC-JC.knee_L).^2),2));
    v_med_lat = cross(v_inf_sup, v_post_ant); % (asse y)
    
    loc_ref.knee_L.x = v_post_ant;
    loc_ref.knee_L.y = v_med_lat;
    loc_ref.knee_L.z = v_inf_sup;
    
    % +++++ Right knee +++++++++++++++++++++++++++++++++
    v_R_THI_KNE = (traj.RTHI-traj.RKNE)./sqrt(sum(((traj.RTHI-traj.RKNE).^2),2));
    v_R_HJC_KNE = (RHJC-traj.RKNE)./sqrt(sum(((RHJC-traj.RKNE).^2),2));
    
    v_post_ant = cross(v_R_HJC_KNE, v_R_THI_KNE)./sqrt(sum(((-v_R_HJC_KNE+v_R_THI_KNE).^2),2));

    JC.knee_R = chordPiG(traj.RTHI,RHJC,traj.RKNE,KneeOS);
    
    % versore che va da RHJC a RKJC (asse z)
    v_inf_sup = (RHJC-JC.knee_R)./sqrt(sum(((RHJC-JC.knee_R).^2),2));
    v_lat_med = cross(v_inf_sup, v_post_ant); % (asse y)

    loc_ref.knee_R.x = v_post_ant;
    loc_ref.knee_R.y = v_lat_med;
    loc_ref.knee_R.z = v_inf_sup;
    
    %% ANKLE
    AnkleOS = (Antro.mDiameter + Antro.ankle_width)/2;
    
    % +++++ Left ankle +++++++++++++++++++++++++++++++++
    v_L_TIB_ANK = (traj.LTIB-traj.LANK)./sqrt(sum(((traj.LTIB-traj.LANK).^2),2));
    v_L_KJC_ANK = (JC.knee_L-traj.LANK)./sqrt(sum(((JC.knee_L-traj.LANK).^2),2));
    v_post_ant = cross(v_L_TIB_ANK, v_L_KJC_ANK)./sqrt(sum(((v_L_KJC_ANK-v_L_TIB_ANK).^2),2));
    
    JC.ankle_L = chordPiG(traj.LTIB,JC.knee_L,traj.LANK,AnkleOS);
    
    % versore che va da LHJC a LKJC (asse z)
    v_inf_sup = (JC.knee_L-JC.ankle_L)./sqrt(sum(((JC.knee_L-JC.ankle_L).^2),2));
    v_med_lat = cross(v_inf_sup, v_post_ant); % (asse y)
    
    loc_ref.ankle_L.x = v_post_ant;
    loc_ref.ankle_L.y = v_med_lat;
    loc_ref.ankle_L.z = v_inf_sup;
    
    % +++++ Right ankle +++++++++++++++++++++++++++++++++
    v_R_TIB_ANK = (traj.RTIB-traj.RANK)./sqrt(sum(((traj.RTIB-traj.RANK).^2),2));
    v_R_KJC_ANK = (JC.knee_R-traj.RANK)./sqrt(sum(((JC.knee_R-traj.RANK).^2),2));
    v_post_ant = cross(v_R_KJC_ANK, v_R_TIB_ANK)./sqrt(sum(((v_R_TIB_ANK-v_R_KJC_ANK).^2),2));
    
    JC.ankle_R = chordPiG(traj.RTIB, JC.knee_R, traj.RANK, AnkleOS);
    
    % versore che va da RHJC a RKJC (asse z)
    v_inf_sup = (JC.knee_R-JC.ankle_R)./sqrt(sum(((JC.knee_R-JC.ankle_R).^2),2));
    v_lat_med = cross(v_inf_sup, v_post_ant);    % (asse y)

    loc_ref.ankle_R.x = v_post_ant;
    loc_ref.ankle_R.y = v_lat_med;
    loc_ref.ankle_R.z = v_inf_sup;

    %% FOOT
    
    % +++++ Left foot +++++++++++++++++++++++++++++++++
    v_L_TOE_HEE = (traj.LHEE-traj.LTOE)./sqrt(sum(((traj.LHEE-traj.LTOE).^2),2)); % (asse z -> blu)
    % vettore di supporto: vettore che collega la caviglia e il ginocchio
    v_L_AJC_KJC = (JC.knee_L-JC.ankle_L)./sqrt(sum(((JC.knee_L-JC.ankle_L).^2),2));

    v_y = cross(v_L_TOE_HEE,v_L_AJC_KJC)./sqrt(sum(((v_L_AJC_KJC-v_L_TOE_HEE).^2),2)); %(asse y -> verde)
    v_x = cross(v_y, v_L_TOE_HEE); %(asse x -> rosso)

    loc_ref.foot_L.x = v_x;
    loc_ref.foot_L.y = v_y;
    loc_ref.foot_L.z = v_L_TOE_HEE;

    % +++++ Right foot +++++++++++++++++++++++++++++++++
    v_R_TOE_HEE = (traj.RHEE-traj.RTOE)./sqrt(sum(((traj.RHEE-traj.RTOE).^2),2)); % (asse z -> blu)
    % vettore di supporto: vettore che collega la caviglia e il ginocchio
    v_R_AJC_KJC = (JC.knee_R-JC.ankle_R)./sqrt(sum(((JC.knee_R-JC.ankle_R).^2),2));

    v_y = cross(v_R_TOE_HEE, v_R_AJC_KJC)./sqrt(sum(((v_R_AJC_KJC-v_R_TOE_HEE).^2),2)); %(asse y -> verde)
    v_x = cross(v_y, v_R_TOE_HEE); %(asse x -> rosso)

    loc_ref.foot_R.x = v_x;
    loc_ref.foot_R.y = v_y;
    loc_ref.foot_R.z = v_R_TOE_HEE;
    
end


