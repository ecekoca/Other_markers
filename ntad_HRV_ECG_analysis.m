%% Heart Rate Variability (HRV) analysis of NTAD using ECG during resting state -eyes open
% (Ece Kocagoncu 2020)

clear, clc
close all

addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/;
addpath /imaging/ek01/ntad/scripts/HRV_scripts/MarcusVollmer-HRV-58badf9/;
ana_dir = '/imaging/projects/cbu/ntad/meg_data/';
task='rest_open';
outdir = '/imaging/ek01/ntad/HRV/';
session = 'BL';
basefile ='transdef';

subjects={'C1001','C1002','C1003','C1004','C1006','C1007','C1008',...
    'C1010','C1011','C1012','C1013','C1014','C1016','C1017',...
    'P1001','P1002','P1004','P1005','P1007','P1008','P1009','P1010',...
    'P1011','P1012','P1015','P1016','P1020','P1021','P1022','P1023',...
    'P1024','P1026','P1027','P1029','P1030','P1031','P1032','P1035',...
    'P1036','P1037','P1038','P1042','P1043','P1045','P1048','P1049',...
    'P1053','P1054','P1055','P1056','P1058','P1060','P1062','P1063',...
    'P1064','P1065','P1066','P1067','P1070','P3002','P3004','P3005','P3006'};%REST_OPEN

% subjects={'C1001','C1002','C1003','C1004','C1006','C1007','C1008',...
%     'C1010','C1011','C1012','C1013','C1014','C1016','C1017',...
%     'P1001','P1002','P1004','P1005','P1007','P1008','P1009','P1010',...
%     'P1011','P1012','P1015','P1016','P1020','P1021','P1022','P1023',...
%     'P1026','P1027','P1029','P1030','P1031','P1032','P1035',...
%     'P1036','P1037','P1038','P1042','P1043','P1045','P1048','P1049',...
%     'P1053','P1054','P1055','P1056','P1058','P1060','P1062','P1063',...
%     'P1064','P1065','P1066','P1067','P1070','P3002','P3004','P3006'};%REST_CLOSED


% C1015 does not have ECG
run_filter=0; % run artefact rejection and bandpass filter

%acer=[96;92;86;92;94;98;85;83;98;92;92;95;97;92;99;77;80;75;72;89;61;84;82;68;93;72;70;85;88;56;64;83;50;50;95;75;87;63;77;81;85;47;83;82;85;96;82;42;75;60;80;77;71;75;61;87;54;57;70;62];
%mmse=[30;28;30;28;30;30;29;30;29;30;29;29;29;30;29;25;27;26;25;29;21;27;27;20;30;23;23;28;28;21;23;27;23;19;30;23;29;20;29;28;26;18;29;25;28;28;27;17;27;22;27;23;24;24;25;27;18;20;22;24];

%%

% Prepare data and save separately
if ~exist([outdir '/ECG/'])
    mkdir([outdir '/ECG/'])
end

if run_filter
    opt.modalities={'ECG'};
    parfor s = 1:length(subjects) % this bit is slow, better to run only once.
        
        ses_dir = [ana_dir subjects{s} '/' session '/' task '/'];
        D = spm_eeg_load([ses_dir 'transdef.mat']);
        
        % bandpass filter the data to the QRS spectrum
        D2=osl_filter(D,[5 15]);
       
        % remove artefacts from the ECG channel
        D3=osl_detect_artefacts(D2,'badchannels',false,opt); %marks bad segments, doesn't remove them!
        ECG=D3(D3.indchantype('ECG'),:,:);
        %ECG_all{s,1}=ECG(good_samples(D3)); % get good segments only:
        %using good segments artficially increases the HRV. 
        ECG_all{s,1}=ECG; % 
    end
    for s=1:length(subjects)
        ECG=ECG_all{s,1};
        save([outdir 'ECG/' subjects{s} '_' session '_' task '_ECG.mat'],'ECG');
    end
end

cd([outdir 'ECG/']);
for s = 1:length(subjects)
    
    load([subjects{s} '_' session '_' task '_ECG.mat']);
    
    % smooth
    d=smooth(ECG,'moving',500);
    
    % remove data points that don't fit the mean
    %d(find(isoutlier(d,'quartile')))=[];
    d(find(isoutlier(d,'mean')))=[];
    
    % baseline correct every second, to correct for trends
    for i=1:1000:(length(d)-1000)
        d(i:i+1000)=d(i:i+1000)-nanmean(d(i:i+1000));
    end; d=d(1:length(d)-1000);
    
    % find the flat parts of the signal and cut them
    aa = diff(d);
    inds = find(abs(aa)< 0.0001);
    d(inds)= [];
    d=-d;
    
    % normalise
    d=mat2gray(d);
    
    % flip if direction is incorrect
    [~,loc1,~,~]=findpeaks(d,'MinPeakDistance',40,'MinPeakHeight',0.8);
    [~,loc2,~,~]=findpeaks(-d,'MinPeakDistance',40,'MinPeakHeight',-0.2);
    if length(loc2)<length(loc1); d=mat2gray(-d); end
    
    % make the rest of the data flat. only peaks remain
    %d(find(d<0.7))=0.1;
    
    % detect R waves
    [~,loc,~,~]=findpeaks(d,'MinPeakDistance',500,'MinPeakHeight',0.8);
    
    %figure('Renderer','painters','Position',[10 10 1500 300]); plot(d(50000:60000));
    %findpeaks(d(1:10000),'MinPeakDistance',500);title(subjects{s})
    %print(gcf,'-dbmp', '-r0',[outdir 'ECG/' subjects{s} '_peaks.bmp'] ,'-r200'); close(gcf)
    
    %average heart beat in the whole recording ~5 mins
    num_pks = numel(loc);
    beats_per_min(s,1) = (num_pks*60)/(length(d)/1000);
    
    % HRV measures
    RMSSD(s,1)= rms(diff(diff(loc)));
    lnRMSSD(s,1)=log(rms(diff(diff(loc)))); % in general lower than 10 (Nunan et al 2010)
    SDNN(s,1)=std(diff(loc));
    
end

tbl=table(subjects',beats_per_min,RMSSD,lnRMSSD,SDNN,'VariableNames',{'SUB','BPM','RMSSD','lnRMSSD','SDNN'}); disp(tbl) % show all the results in a table
save([outdir '/HRV_metrics_baseline_rest_open.mat'],'tbl'); % save
%[r,p]=corrcoef(RMSSD,beats_per_min)
%figure; scatter(RMSSD,beats_per_min);lsline



