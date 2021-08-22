%% Plasma biomarker analyses (Ece K 2021)

clear all; close all; clc
outdir='/imaging/rowe/users/ek01/prv/ntad/plasma_ptau181/';
fname='input.mat';

yellow=[255 193 0]./255;
yellow2=[245 222 159]./255;
grey=[170 170 170]./255;
grey2=[200 200 200]./255;

%%

cd(outdir)
load(fname);
subs=T.SUBS;

c_ind=intersect(find(contains(subs,'C')),find(T.AMY_ST==0));
p_ind=intersect(find(contains(subs,'P')),find(T.AMY_ST==1));

ptau=T.PTAU181; amy=T.TAU_AB; amy2=T.FBB_SUVR; 
age=T.AGE; mmse=T.MMSE; acer=T.ACER; cdr=T.CDR; sob=T.CDR_SOB;
acer_mem=T.ACER_MEM;hip=T.HIP; gmv=T.GMV; wmv=T.WMV; edu=T.YOE; sex=T.SEX;
tem_pet=T.TEMP_PET; hip_pet=T.HIP_PET; pl_ab40=T.AB40; pl_ab42=T.AB42; pl_ab_rat=T.AB42_40;
pl_gfap=T.GFAP; pl_nfl=T.NFL; pl_ptau_ab=T.PTAU181_AB42;

%% Boxplot

m_col=[0 192 163]./255;
a_col=[63 14 137]./255;
c_col=[255 225 0]./255;
var=[ptau pl_ab40 pl_ab42 pl_ab_rat pl_gfap pl_nfl pl_ptau_ab];
varnm={'P-TAU181','AB40','AB42','AB42/40','GFAP','NFL','PTAU181/AB42'};

groups = [0.1*ones(size(c_ind')), 0.4*ones(size(p_ind'))]';

for v=1:size(var,2)
    figure('color','w'); set(gca,'FontSize',14,'Color','w')
    H=notBoxPlot([var(c_ind,v); var(p_ind,v)],groups,'jitter',0.2)
    set([H.data],'MarkerSize',4,...
        'markerFaceColor',[1 1 1]*0.25,...
        'markerEdgeColor', 'none')
    set([H(1:2).semPtch],'FaceColor',grey,'EdgeColor','none')
    set([H(1:2).sdPtch],'FaceColor',grey2,'EdgeColor','none')
    set([H.mu],'Color','k');
    set([H(1).data],'MarkerSize',7,'markerFaceColor',c_col,'markerEdgeColor', 'none')
    set([H(2).data],'MarkerSize',7,'markerFaceColor',a_col,'markerEdgeColor', 'none')
    set(gca, 'XtickLabel', {'C','MCI/AD'}); %ylim([0 4.5]);
    ylabel([varnm{v} ' (pg/ml)']); title ('Group differences');
    print(gcf,[outdir 'plasma_' strrep(varnm{v},'/','_') '_boxplots.bmp'],'-dbmp','-r300'); close(gcf)
end


mci_ind=intersect(p_ind,find(strcmp(T.GROUP,'MCI')));
ad_ind=intersect(p_ind,find(strcmp(T.GROUP,'AD')));

groups = [0.1*ones(size(c_ind')), 0.4*ones(size(mci_ind')),  0.7*ones(size(ad_ind'))]';

for v=1:size(var,2)
    figure('color','w'); set(gca,'FontSize',14,'Color','w')
    H=notBoxPlot([var(c_ind,v); var(mci_ind,v); var(ad_ind,v)],groups,'jitter',0.2) %
    set([H([1:3]).semPtch],'FaceColor',grey,'EdgeColor','none')
    set([H([1:3]).sdPtch],'FaceColor',grey2,'EdgeColor','none')
    set([H.mu],'Color','k');
    set([H(1).data],'MarkerSize',7,'markerFaceColor',c_col,'markerEdgeColor', 'none')
    set([H(2).data],'MarkerSize',7,'markerFaceColor',m_col,'markerEdgeColor', 'none')
    set([H(3).data],'MarkerSize',7,'markerFaceColor',a_col,'markerEdgeColor', 'none')
    set(gca, 'XtickLabel', {'C','MCI', 'AD'}); %ylim([0 4.5]);
    ylabel([varnm{v} ' (pg/ml)']); title ('Group differences');
    print(gcf,[outdir 'plasma_' strrep(varnm{v},'/','_') '_boxplots_3grp.bmp'],'-dbmp','-r300'); close(gcf)
end

p1_ind=intersect(p_ind,find(T.MMSE<25));
p2_ind=intersect(p_ind,find(T.MMSE>25));

groups = [0.1*ones(size(c_ind')), 0.4*ones(size(p1_ind')),  0.7*ones(size(p2_ind'))]';

for v=1:size(var,2)
    figure('color','w'); set(gca,'FontSize',14,'Color','w')
    H=notBoxPlot([var(c_ind,v); var(p2_ind,v); var(p1_ind,v)],groups,'jitter',0.2) %
    set([H([1:3]).semPtch],'FaceColor',grey,'EdgeColor','none')
    set([H([1:3]).sdPtch],'FaceColor',grey2,'EdgeColor','none')
    set([H.mu],'Color','k');
    set([H(1).data],'MarkerSize',7,'markerFaceColor',c_col,'markerEdgeColor', 'none')
    set([H(2).data],'MarkerSize',7,'markerFaceColor',m_col,'markerEdgeColor', 'none')
    set([H(3).data],'MarkerSize',7,'markerFaceColor',a_col,'markerEdgeColor', 'none')
    set(gca, 'XtickLabel', {'C','>24', '<24'}); %ylim([0 4.5]);
    ylabel([varnm{v} ' (pg/ml)']); title ('Group differences');
    print(gcf,[outdir 'plasma_' strrep(varnm{v},'/','_') '_boxplots_3grp_by_MMSE.bmp'],'-dbmp','-r300'); close(gcf)
end

var=[T.P300 T.P300_SLOPE T.MMN_LAT T.P300_QSLOPE T.EO_FREQ T.EC_FREQ T.EO_POW T.EC_POW];
varnm={'P300','P300 Slope','MMN Latency','P300 Q Slope','Eyes open Alpha Fr',...
    'Eyes closed Alpha Fr','Eyes open Alpha Pow','Eyes closed Alpha Pow'};

groups = [0.1*ones(size(c_ind')), 0.4*ones(size(p_ind'))]';

for v=1:size(var,2)
    figure('color','w'); set(gca,'FontSize',14,'Color','w')
    H=notBoxPlot([var(c_ind,v); var(p_ind,v)],groups,'jitter',0.2)
    set([H.data],'MarkerSize',4,...
        'markerFaceColor',[1 1 1]*0.25,...
        'markerEdgeColor', 'none')
    set([H(1:2).semPtch],'FaceColor',grey,'EdgeColor','none')
    set([H(1:2).sdPtch],'FaceColor',grey2,'EdgeColor','none')
    set([H.mu],'Color','k');
    set([H(1).data],'MarkerSize',7,'markerFaceColor',c_col,'markerEdgeColor', 'none')
    set([H(2).data],'MarkerSize',7,'markerFaceColor',a_col,'markerEdgeColor', 'none')
    set(gca, 'XtickLabel', {'C','MCI/AD'}); %ylim([0 4.5]);
    ylabel([varnm{v}]); title ('Group differences');
    print(gcf,[outdir 'MEG_' strrep(varnm{v},'/','_') '_boxplots.bmp'],'-dbmp','-r300'); close(gcf)
end

cn=length(c_ind); pn=length(p_ind); mcin=length(mci_ind); adn=length(ad_ind);

for v=1:size(var,2)
    X=[]; X(1:cn,1)=1; X(cn+1:cn+pn,2)=1; %X(:,3)=zscore(age([c_ind; p_ind]));
    y=var([c_ind; p_ind],v); idx=find(isnan(y));y(idx)=[]; X(idx,:)=[];
    y=zscore(y);
    [t,F,p,df,R2,cR2,B] = glm_rik(y,X,[1 -1 ]',1);
    results(v,1)=t; results(v,2)=p;
end
%% Group differences
cn=length(c_ind); pn=length(p_ind); mcin=length(mci_ind); adn=length(ad_ind);

for v=1:size(var,2)
    X=[]; X(1:cn,1)=1; X(cn+1:cn+pn,2)=1; X(:,3)=zscore(age([c_ind; p_ind]));
    y=zscore(var([c_ind; p_ind],v));
    [t,F,p,df,R2,cR2,B] = glm_rik(y,X,[1 -1 0]',1);
    results(v,1)=t; results(v,2)=p;
end

for v=1:size(var,2)
    X=[]; X(1:cn,1)=1; X(cn+1:cn+mcin,2)=1; X(:,3)=zscore(age([c_ind; mci_ind]));
    y=zscore(var([c_ind; mci_ind],v));
    [t,F,p,df,R2,cR2,B] = glm_rik(y,X,[1 -1 0]',1);
    results(v,3)=t; results(v,4)=p;
end

for v=1:size(var,2)
    X=[]; X(1:cn,1)=1; X(cn+1:cn+adn,2)=1; X(:,3)=zscore(age([c_ind; ad_ind]));
    y=zscore(var([c_ind; ad_ind],v));
    [t,F,p,df,R2,cR2,B] = glm_rik(y,X,[1 -1 0]',1);
    results(v,5)=t; results(v,6)=p;
end

for v=1:size(var,2)
    X=[]; X(1:mcin,1)=1; X(mcin+1:mcin+adn,2)=1; X(:,3)=zscore(age([mci_ind; ad_ind]));
    y=zscore(var([mci_ind; ad_ind],v));
    [t,F,p,df,R2,cR2,B] = glm_rik(y,X,[1 -1 0]',1);
    results(v,7)=t; results(v,8)=p;
end

%% ROC

for v=1:size(var,2)
    
    d1=roc_curve(var(c_ind,v),var(p_ind,v));
    lines{1}(:,:)=d1.curve; auc(1)=d1.param.AROC;
    d2=roc_curve(var(c_ind,v),var(mci_ind,v));
    lines{2}(:,:)=d2.curve; auc(2)=d2.param.AROC;
    d3=roc_curve(var(c_ind,v),var(ad_ind,v));
    lines{3}(:,:)=d3.curve; auc(3)=d3.param.AROC;
    d4=roc_curve(var(mci_ind,v),var(ad_ind,v));
    lines{4}(:,:)=d4.curve; auc(4)=d4.param.AROC;
    for c=1:4; if auc(c)<0.5; auc(c)=1-auc(c); lines{c}(:,:)=1-lines{c}(:,:);end; end
    
    figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.3, 0.4])
    plot(1-lines{3}(:,2),lines{3}(:,1),'Color',[255 225 0]./255 ,'LineWidth',2.5); hold on;
    plot(1-lines{2}(:,2),lines{2}(:,1),'Color',[0.6 0.6 0.6],'LineWidth',2.5);
    plot(1-lines{1}(:,2),lines{1}(:,1),'Color',[63 14 137]./255,'LineWidth',2.5);
    plot(1-lines{4}(:,2),lines{4}(:,1),'Color',[0 192 163]./255,'LineWidth',2.5);
    axis square;grid on; xlabel('1 - Specificity'); ylabel('Sensitivity');
    legend(['C x AD: ' num2str(auc(3))],['C x MCI: ' num2str(auc(2))],...
        ['C x P: ' num2str(auc(1))],['MCI x AD: ' num2str(auc(4))],'Location','SouthEast')
    set(gca,'FontSize',14); title(varnm{v});
    print(gcf,[outdir 'ROC_curves_plasma_' strrep(varnm{v},'/','_')  '.bmp'],'-dbmp','-r300'); close(gcf)
    
end; clear lines

for v=1:size(var,2)
    d1=roc_curve(var(c_ind,v),var(p_ind,v));
    lines(v,:,:)=d1.curve; auc(v)=d1.param.AROC; 
end
for c=1:size(var,2); if auc(c)<0.5; auc(c)=1-auc(c); lines(c,:,:)=1-lines(c,:,:);end; end
figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.3, 0.4]); hold on;
for v=1:size(var,2); plot(1-lines(v,:,2),lines(v,:,1),'LineWidth',2.5); end
axis square;grid on; xlabel('1 - Specificity'); ylabel('Sensitivity');
legend([varnm{1} ': ' num2str(auc(1))],[varnm{2} ': ' num2str(auc(2))],...
    [varnm{3} ': ' num2str(auc(3))],[varnm{4} ': ' num2str(auc(4))],...
    [varnm{5} ': ' num2str(auc(5))],[varnm{6} ': ' num2str(auc(6))],...
    [varnm{7} ': ' num2str(auc(7))],'Location','EastOutside')
set(gca,'FontSize',14); title('Plasma markers C x P');
print(gcf,[outdir 'ROC_curves_plasma_all_markers_C-P.bmp'],'-dbmp','-r300'); close(gcf)

clear lines auc
for v=1:size(var,2)
    d1=roc_curve(var(mci_ind,v),var(ad_ind,v));
    lines(v,:,:)=d1.curve; auc(v)=d1.param.AROC; 
end
for c=1:size(var,2); if auc(c)<0.5; auc(c)=1-auc(c); lines(c,:,:)=1-lines(c,:,:);end; end
figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.3, 0.4]); hold on;
for v=1:size(var,2); plot(1-lines(v,:,2),lines(v,:,1),'LineWidth',2.5); end
axis square;grid on; xlabel('1 - Specificity'); ylabel('Sensitivity');
legend([varnm{1} ': ' num2str(auc(1))],[varnm{2} ': ' num2str(auc(2))],...
    [varnm{3} ': ' num2str(auc(3))],[varnm{4} ': ' num2str(auc(4))],...
    [varnm{5} ': ' num2str(auc(5))],[varnm{6} ': ' num2str(auc(6))],...
    [varnm{7} ': ' num2str(auc(7))],'Location','EastOutside')
set(gca,'FontSize',14); title('Plasma markers MCI x AD');
print(gcf,[outdir 'ROC_curves_plasma_all_markers_MCI-AD.bmp'],'-dbmp','-r300'); close(gcf)

%% SVM (and ROC)

var=[T.MMSE T.ACER T.ACER_MEM T.CDR_SOB ...
    T.PTAU181 T.AB40 T.AB42 T.AB42_40 T.GFAP T.NFL T.PTAU181_AB42 ...
    T.TTAU T.AB1_42 T.TAU_AB T.PTAU T.TTAU_PTAU ...
    T.FBB_SUVR T.TEMP_PET T.HIP_PET ...
    T.GMV T.WMV T.HIP T.HIP_CA1 T.HIP_SUB ...
    T.MMN_LAT T.P300 T.P300_SLOPE T.P300_QSLOPE T.EO_FREQ T.EC_FREQ T.EO_POW T.EC_POW];
varnm={'MMSE','ACER','ACER memory','CDR SOB',...
    'Plasma ptau181','Plasma AB40','Plasma AB42','Plasma AB42/40', 'Plasma GFAP','Plasma NFL','Plasma ptau181/AB42',...
    'CSF total tau','CSF AB1-42','CSF Tau/AB1-42','CSF ptau','CSF ttau/ptau',...
    'Global SUVR','Temporal SUVR','Hip SUVR',...
    'GMV','WMV','Hip GMV','Hip CA1 GM','Hip sub GM',...
    'MMN latency','P300 mean','P300 slope','P300 q slope','Eyes open freq','Eyes closed freq',...
    'Eyes open power','Eyes closed power'};

% Controls vs patients plasma 
cl(1:length(c_ind),1)=1;cl(length(c_ind)+1:length(c_ind)+length(p_ind),1)=2;
for v=1:size(var,2)
    A=[var([c_ind; p_ind],v) cl];
    idx=find(isnan(A(:,1))); A(idx,:)=[];
    try
        acc=repeated_CV(A,5,250);
        m_acc(v,1)=mean(acc);
    catch
        m_acc(v,1)=NaN;
    end
    d1=roc_curve(A(find(A(:,2)==1),1),A(find(A(:,2)==2),1));
    k=d1.param.AROC; if k<0.5; k=1-k; end; auc(v,1)=k;
end

% MCI vs AD plasma
clear A cl
cl(1:length(mci_ind),1)=1;cl(length(mci_ind)+1:length(mci_ind)+length(ad_ind),1)=2;
for v=1:size(var,2)
    A=[var([mci_ind; ad_ind],v) cl]; 
    idx=find(isnan(A(:,1))); A(idx,:)=[]; 
    try
        acc=repeated_CV(A,5,250);
        m_acc2(v,1)=mean(acc);
    catch
        m_acc2(v,1)=NaN;
    end    
    d1=roc_curve(A(find(A(:,2)==1),1),A(find(A(:,2)==2),1));
    k=d1.param.AROC; if k<0.5; k=1-k; end; auc2(v,1)=k;
end

T2=table(varnm',auc,m_acc,auc2,m_acc2,'VariableNames',{'METRICS','AUC1','ACC1','AUC2','ACC2'})
save([outdir 'SVM_AUC_results.mat'],'T2')

%% Plot correlation matrix

var=[T.SEX T.AGE T.YOE T.MMSE T.ACER T.ACER_MEM T.CDR_SOB T.PTAU181 T.AB42_40 T.GFAP T.NFL T.PTAU181_AB42...
    T.TTAU T.AB1_42 T.TAU_AB T.PTAU T.TTAU_PTAU T.FBB_SUVR T.TEMP_PET T.HIP_PET ...
    T.GMV T.WMV T.HIP T.HIP_CA1 T.HIP_SUB ...
    T.MMN_LAT T.P300 T.P300_SLOPE T.P300_QSLOPE T.EO_FREQ T.EC_FREQ T.EO_POW T.EC_POW];
varnm={'Sex','Age','YoE','MMSE','ACER','ACER memory','CDR SOB','Plasma ptau181','Plasma AB42/40',...
    'Plasma GFAP','Plasma NFL','Plasma ptau181/AB42','CSF total tau','CSF AB1-42','CSF Tau/AB1-42','CSF ptau','CSF ttau/ptau',...
    'Global SUVR','Temporal SUVR','Hip SUVR','GMV','WMV','Hip GMV','Hip CA1 GM','Hip sub GM',...
    'MMN latency','P300 mean','P300 slope','P300 q slope','Eyes open freq','Eyes closed freq',...
    'Eyes open power','Eyes closed power'};

for i=1:size(var,2)
    for j=1:size(var,2)
        [r p]=corrcoef(var([c_ind; p_ind],i),var([c_ind; p_ind],j),'rows','complete');
        rmat(i,j)=r(1,2); pmat(i,j)=p(1,2)/2; %one tailed
        [r p]=corrcoef(var([p_ind],i),var([p_ind],j),'rows','complete');
        rmat2(i,j)=r(1,2); pmat2(i,j)=p(1,2)/2; %one tailed
        [r p]=corrcoef(var([c_ind],i),var([c_ind],j),'rows','complete');
        rmat3(i,j)=r(1,2); pmat3(i,j)=p(1,2)/2; %one tailed
    end
end
rmat(find(isnan(rmat)))=0; rmat2(find(isnan(rmat2)))=0; rmat3(find(isnan(rmat3)))=0;
figure('color','w'); imagesc(rmat);
axis square; colormap('redblue'); caxis([-0.5 0.5]); colorbar
xticks(1:size(var,2));yticks(1:size(var,2))
set(gca,'XTickLabel',varnm,'YTickLabel',varnm); xtickangle(40)
[a,b]=ind2sub(size(pmat),find(pmat<=0.05)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',4,'markeredgecolor','w','markerfacecolor','w')
end
[a,b]=ind2sub(size(pmat),find(pmat<=0.01)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',6,'markeredgecolor','w','markerfacecolor','w')
end; title('Marker intercorrelations')
print(gcf,'-dbmp', '-r0',[outdir 'correlation_matrix_all.bmp'],'-r200');% close(gcf)

% Plot plasma correlations only
plasma_rmat=rmat(8:12,:); plasma_pmat=pmat(8:12,:);
figure('color','w'); imagesc(plasma_rmat);
axis equal; colormap('redblue'); caxis([-0.5 0.5]); colorbar
xticks(1:size(var,2));yticks(1:5)
set(gca,'XTickLabel',varnm,'YTickLabel',varnm(8:12)); xtickangle(40)
[a,b]=ind2sub(size(plasma_pmat),find(plasma_pmat<=0.05)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',4,'markeredgecolor','w','markerfacecolor','w')
end
[a,b]=ind2sub(size(plasma_pmat),find(plasma_pmat<=0.01)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',6,'markeredgecolor','w','markerfacecolor','w')
end; ylim([0.5 5.5])
print(gcf,'-dbmp', '-r0',[outdir 'correlation_matrix_plasma.bmp'],'-r200');% close(gcf)

plasma_rmat=rmat2(8:12,:); plasma_pmat=pmat2(8:12,:);
figure('color','w'); imagesc(plasma_rmat);
axis equal; colormap('redblue'); caxis([-0.5 0.5]); colorbar
xticks(1:size(var,2));yticks(1:5)
set(gca,'XTickLabel',varnm,'YTickLabel',varnm(8:12)); xtickangle(40)
[a,b]=ind2sub(size(plasma_pmat),find(plasma_pmat<=0.05)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',4,'markeredgecolor','w','markerfacecolor','w')
end
[a,b]=ind2sub(size(plasma_pmat),find(plasma_pmat<=0.01)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',6,'markeredgecolor','w','markerfacecolor','w')
end; ylim([0.5 5.5])
print(gcf,'-dbmp', '-r0',[outdir 'correlation_matrix_plasma_patients_only.bmp'],'-r200');% close(gcf)

mat=[rmat(8,:);rmat3(8,:);rmat2(8,:)]; %ptau181 row
mat(find(isnan(mat)))=0;
mat2=[pmat(8,:);pmat3(8,:);pmat2(8,:)];
figure('Color','w'); imagesc(mat)
axis equal; colormap('redblue'); caxis([-0.5 0.5]); colorbar
xticks(1:size(var,2));yticks(1:3)
set(gca,'XTickLabel',varnm,'YTickLabel',{'Controls & Patients','Controls','Patients'}); xtickangle(40)
[a,b]=ind2sub(size(mat2),find(mat2<=0.05)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',4,'markeredgecolor','w','markerfacecolor','w')
end
[a,b]=ind2sub(size(mat2),find(mat2<=0.01)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',6,'markeredgecolor','w','markerfacecolor','w')
end; ylim([0.5 3.5]); title('Ptau181 correlations')
print(gcf,'-dbmp', '-r0',[outdir 'correlation_matrix_ptau181_diff_groups.bmp'],'-r200');% close(gcf)


% clustergram(rmat,'Colormap','redblue','RowLabels',varnm,'ColumnLabels',varnm,...
%     'ColumnLabelsRotate',40,'DisplayRange',0.5); 
% print(gcf,'-dbmp', '-r0',[outdir 'clustergram_all.bmp'],'-r200');% close(gcf)
% 
% tree=linkage(rmat,'average');
% h=dendrogram(tree,0,'ColorThreshold',2.5,'labels',varnm);
% set(h,'LineWidth',1.4); xtickangle(40); box on
% print(gcf,'-dbmp', '-r0',[outdir 'dendrogram_all.bmp'],'-r200');% close(gcf)

% Plot MEG correlations only
meg_rmat=rmat(26:33,:); meg_pmat=pmat(26:33,:);
figure('color','w'); imagesc(meg_rmat);
axis equal; colormap('redblue'); caxis([-0.5 0.5]); colorbar
xticks(1:size(var,2));yticks(1:8)
set(gca,'XTickLabel',varnm,'YTickLabel',varnm(26:33)); xtickangle(40)
[a,b]=ind2sub(size(meg_pmat),find(meg_pmat<=0.05)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',4,'markeredgecolor','w','markerfacecolor','w')
end
[a,b]=ind2sub(size(meg_pmat),find(meg_pmat<=0.01)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',6,'markeredgecolor','w','markerfacecolor','w')
end; ylim([0.5 8.5])
print(gcf,'-dbmp', '-r0',[outdir 'correlation_matrix_meg.bmp'],'-r200');% close(gcf)

meg_rmat=rmat2(26:33,:); meg_pmat=pmat2(26:33,:);
figure('color','w'); imagesc(meg_rmat);
axis equal; colormap('redblue'); caxis([-0.5 0.5]); colorbar
xticks(1:size(var,2));yticks(1:8)
set(gca,'XTickLabel',varnm,'YTickLabel',varnm(26:33)); xtickangle(40)
[a,b]=ind2sub(size(meg_pmat),find(meg_pmat<=0.05)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',4,'markeredgecolor','w','markerfacecolor','w')
end
[a,b]=ind2sub(size(meg_pmat),find(meg_pmat<=0.01)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',6,'markeredgecolor','w','markerfacecolor','w')
end; ylim([0.5 8.5])
print(gcf,'-dbmp', '-r0',[outdir 'correlation_matrix_meg_patients_only.bmp'],'-r200');% close(gcf)


% Partial correlations age corrected
for i=1:size(var,2)
    for j=1:size(var,2)
        try
            [r p]=partialcorr(var([c_ind; p_ind],i),var([c_ind; p_ind],j),T.AGE([c_ind;p_ind]),'rows','complete');
            rmat(i,j)=r; pmat(i,j)=p;
        catch
            rmat(i,j)=0; pmat(i,j)=1;
        end
    end
end
rmat(find(isnan(rmat)))=0;
figure('color','w'); imagesc(rmat);
axis square; colormap('redblue'); caxis([-0.5 0.5]); colorbar
xticks(1:size(var,2));yticks(1:size(var,2))
set(gca,'XTickLabel',varnm,'YTickLabel',varnm); xtickangle(40)
[a,b]=ind2sub(size(pmat),find(pmat<=0.05)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',4,'markeredgecolor','w','markerfacecolor','w')
end
[a,b]=ind2sub(size(pmat),find(pmat<=0.01)); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',6,'markeredgecolor','w','markerfacecolor','w')
end; title('Marker intercorrelations age cor')
print(gcf,'-dbmp', '-r0',[outdir 'correlation_matrix_age_cor.bmp'],'-r200');% close(gcf)

tree=linkage(rmat,'average');
h=dendrogram(tree,0,'ColorThreshold',2.5,'labels',varnm);
set(h,'LineWidth',1.4); xtickangle(40); box on
print(gcf,'-dbmp', '-r0',[outdir 'dendrogram_all_age_cor.bmp'],'-r200');% close(gcf)

%% Plot scatterplots

c_col=[255 225 0]./255;
m_col=[0 192 163]./255;
a_col=[63 14 137]./255;

pvar=[T.PTAU181 T.AB40 T.AB42 T.AB42_40 T.GFAP T.NFL T.PTAU181_AB42];
var=[T.AGE T.YOE T.MMSE T.ACER T.ACER_MEM T.CDR_SOB ...
    T.TTAU T.AB1_42 T.TAU_AB T.PTAU T.TTAU_PTAU ...
    T.FBB_SUVR T.TEMP_PET T.HIP_PET ...
    T.GMV T.WMV T.HIP T.HIP_CA1 T.HIP_SUB ...
    T.MMN_LAT T.P300 T.P300_SLOPE T.P300_QSLOPE T.EO_FREQ T.EC_FREQ T.EO_POW T.EC_POW];
pvarnm={'Plasma ptau181','Plasma AB40','Plasma AB42','Plasma AB42/40', 'Plasma GFAP','Plasma NFL','Plasma ptau181/AB42'};
varnm={'Age','YoE','MMSE','ACER','ACER memory','CDR SOB',...
    'CSF total tau','CSF AB1-42','CSF Tau/AB1-42','CSF ptau','CSF ttau/ptau',...
    'Global SUVR','Temporal SUVR','Hip SUVR',...
    'GMV','WMV','Hip GMV','Hip CA1 GM','Hip sub GM',...
    'MMN latency','P300 mean','P300 slope','P300 q slope','Eyes open freq','Eyes closed freq',...
    'Eyes open power','Eyes closed power'};

for i=1:size(pvar,2)
    for j=1:size(var,2)
        d1=pvar(:,i);d2=var(:,j);
        idx=unique([find(isoutlier(pvar(:,i))); find(isoutlier(var(:,j)))]);
        d1(idx)=NaN; d2(idx)=NaN;
        figure('color','w');
        scatter(d1([c_ind;mci_ind;ad_ind]),d2([c_ind;mci_ind;ad_ind]),80,'w','filled');
        h=lsline; h.LineWidth=3; h.Color='k';hold on;
        scatter(d1(c_ind),d2(c_ind),120,c_col,'filled');
        scatter(d1(mci_ind),d2(mci_ind),120,m_col,'filled');
        scatter(d1(ad_ind),d2(ad_ind),120,a_col,'filled'); axis square; box on;
        xlabel([pvarnm{i} ' (pg/ml)']); ylabel(varnm{j});
        set(gca,'FontSize',14); grid on;
        print(gcf,[outdir 'scatterplot_' strrep(strrep(pvarnm{i},' ','_'),'/','_') '_' strrep(strrep(varnm{j},' ','_'),'/','_') '.bmp'],'-dbmp','-r300'); close(gcf)
    end
end


