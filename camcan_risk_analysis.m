%% Test risk factors in CamCAN3000 using chi square (Ece K, 2020)

close all
clear all; clc

addpath('/imaging/ek01/camcan3000_EK/scripts/')
addpath('/imaging/ek01/camcan3000_EK/')
outdir='/imaging/ek01/camcan3000_EK/risk_factors/';
fname='input.mat';
c_color=[0 84 147]./255; f_color=[193 193 193]./255; %[255 65 81]./255

%%

cd(outdir)
load(fname)

prop_sum=T.C+T.F; ind=find(prop_sum==0); T(ind,:)=[];
varnm=T.risks;
risk_type=T.risk_type;
props=[T.C T.F]; cs=T.CN; fs=T.FN; clear T

for v=1:length(varnm)
   [h,p,chistat,df]=prop_test(props(v,:),[cs(v) fs(v)],false);
   chi2(v,1)=chistat; pval_2tail(v,1)=p; pval_1tail(v,1)=p/2;
end

ratio_c=props(:,1)./cs*100; %percentage of controls where event is true
ratio_f=props(:,2)./fs*100; %percentage of frail where event is true

c_f_OR=(props(:,1)./props(:,2))./((cs-props(:,1))./(fs-props(:,2))); %Control/Frail odds ratio
f_c_OR=1./c_f_OR; %Frail/Control odds ratio

for v=1:length(varnm)
   se(v,1)=1.96*(sqrt((1/props(v,1))+(1/props(v,2))+(1/(cs(v)-props(v,1)))+(1/(fs(v)-props(v,2)))));
   cf_OR_CI(v,:)=[exp(log(c_f_OR(v))-se(v,1)) exp(log(c_f_OR(v))+se(v,1))];
   fc_OR_CI(v,:)=[exp(log(f_c_OR(v))-se(v,1)) exp(log(f_c_OR(v))+se(v,1))];
end


T1=table(varnm,ratio_c,ratio_f,c_f_OR,f_c_OR,chi2,pval_2tail,pval_1tail,'VariableNames',{'Risk_factors','Control_perc','Frail_perc','C_F_OR','F_C_OR','Chi2','p_2tail','p_1tail'})

save([outdir 'risk_factors_chi2_results.mat'],'T1');

%% Plot bar plots

T2=T1;
ind=find(T2.p_1tail>0.05); T2(ind,:)=[];
ind=find(T2.F_C_OR>T2.C_F_OR); chi2=T2.Chi2; chi2(ind)=-chi2(ind);
T2.Chi2=chi2;
T2=sortrows(T2,6,'ascend'); x=1:length(T2.Risk_factors);

c=customcolormap([0 0.5 1],[0 0 1; 0.8 0.8 0.8; 1 0 0]); k1=length(find(chi2<0));
st=floor(128/k1); cmap=c(1:st:129,:);
k2=length(find(chi2>0)); st=floor(128/k2);
cmap2=c(129:st:end,:);cmap=[cmap(1:k1,:); cmap2(1:k2,:)];
OR=[T2.F_C_OR(T2.F_C_OR>T2.C_F_OR); T2.C_F_OR(T2.C_F_OR>T2.F_C_OR)];
OR=round(OR*100)/100; OR(OR==Inf)=1.00; 

figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.3, 0.55])
h=superbar(T2.Chi2,'Orientation','h','BarFaceColor',cmap); box off 
for i=1:length(T2.Risk_factors)
    if T2.p_1tail(i)<0.05 && T2.p_1tail(i)>0.01; pre=''; post=['   (' num2str(T2.p_1tail(i)) ') ']; 
        elseif T2.p_1tail(i)<0.01 && T2.p_1tail(i)>0.001; pre=''; post=['   (' num2str(T2.p_2tail(i)) ') '];
        elseif T2.p_1tail(i)<0.001; pre=''; post=['   (' num2str(T2.p_1tail(i)) ') ']; 
    else pre=''; post=''; end
    if T2.Chi2(i)<0; drc='right';drc2='left'; bgn=2; else; drc='left'; drc2='right';bgn=-2; end
    text(T2.Chi2(i),i,[pre num2str(OR(i)) post],'HorizontalAlignment',drc,'fontsize',9);
    text(bgn,i,[pre T2.Risk_factors{i}],'HorizontalAlignment',drc2,'fontsize',9);
end; 
k=round(max(abs(T2.Chi2)))+5; xlim([-k k]); 
grid on; 
set(gca,'XColor','none'); set(gca,'YTick',[],'YColor','w')
%xtickangle(90)
print(gcf,['/imaging/ek01/camcan3000_EK/all_significant_factors_bar_plots_uncor.bmp'],'-dbmp','-r300');

thresh=0.05/length(T2.Risk_factors);
figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.3, 0.55])
h=superbar(T2.Chi2,'Orientation','h','BarFaceColor',cmap); box off 
for i=1:length(T2.Risk_factors)
    if T2.p_2tail(i)<thresh; pre='\bf'; post=['   (' num2str(T2.p_2tail(i)) ') ']; 
    else pre=''; post=['   (' num2str(T2.p_2tail(i)) ') ']; end
    if T2.Chi2(i)<0; drc='right';drc2='left'; bgn=2; else; drc='left'; drc2='right';bgn=-2; end
    text(T2.Chi2(i),i,[pre num2str(OR(i)) post],'HorizontalAlignment',drc,'fontsize',9);
    text(bgn,i,[pre T2.Risk_factors{i}],'HorizontalAlignment',drc2,'fontsize',9);
end; 
k=round(max(abs(T2.Chi2)))+5; xlim([-k k]); 
grid on; 
set(gca,'XColor','none'); set(gca,'YTick',[],'YColor','w')
%xtickangle(90)
print(gcf,['/imaging/ek01/camcan3000_EK/all_significant_factors_bar_plots_cor.bmp'],'-dbmp','-r300');

types=unique(risk_type);

t=1;
ind=find(strcmp(risk_type,types{t}));
risks=T1.Risk_factors(ind);pval=T1.p_1tail(ind);
C_perc=T1.Control_perc(ind);F_perc=T1.Frail_perc(ind);
CF_OR=T1.C_F_OR(ind);FC_OR=T1.F_C_OR(ind);
chi=T1.Chi2(ind);
sum_perc=C_perc+F_perc;
for i=1:length(risks)
    if C_perc(i)<F_perc(i)
        sum_perc(i)=-sum_perc(i);
        chi(i)=-chi(i);
    end
end
ind=find(sum_perc<0); C_perc(ind)=-C_perc(ind);F_perc(ind)=-F_perc(ind);
T2=table(risks,C_perc,F_perc,CF_OR,FC_OR,sum_perc,chi,pval);
T2=sortrows(T2,7,'ascend'); x=1:length(risks);y=[T2.C_perc T2.F_perc];
c=customcolormap([0 0.5 1],[0 0 1; 0.8 0.8 0.8; 1 0 0]); k1=length(find(chi<0));
st=floor(128/k1); cmap=c(1:st:129,:);
k2=length(find(chi>0)); st=floor(128/k2);
cmap2=c(129:st:end,:);cmap=[cmap(1:k1,:); cmap2(1:k2,:)];

figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.3, 0.6])
h=superbar(T2.chi,'Orientation','v','BarFaceColor',cmap); box off 
endpt=T2.chi; i=find(endpt>0);endpt(i)=endpt(i)+6;
j=find(endpt<0);endpt(j)=endpt(j)-6; 
OR=[T2.FC_OR(T2.FC_OR>T2.CF_OR); T2.CF_OR(T2.CF_OR>T2.FC_OR)];
OR=round(OR*100)/100; OR(OR==Inf)=1.00;OR=abs(OR);
for i=1:length(T2.risks)
    if T2.pval(i)<0.05 && T2.pval(i)>0.01; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')']; 
        elseif T2.pval(i)<0.01 && T2.pval(i)>0.001; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')'];
        elseif T2.pval(i)<0.001; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')']; 
    else pre=''; post=''; end
    
    if endpt(i)<0; drc='right';drc2='left'; bgn=0.5; else; drc='left'; drc2='right';bgn=-0.5; end
    text(i,bgn,[pre T2.risks{i} ' (' num2str(OR(i)) ')'],'HorizontalAlignment',drc2,'Rotation',30,'fontsize',9);
end
k=round(max(abs(endpt))); ylim([-k k]); %legend({'Control','Cognitively Frail'},'NumColumns',2,'Location','SouthOutside')
grid on; %ylabel('Nutrition, smoking alcohol & use','Fontweight','bold'); 
set(gca,'YColor','none'); set(gca,'XTick',[],'XColor','w')
%xtickangle(90)
print(gcf,['/imaging/ek01/camcan3000_EK/nutrition_bar_plots.bmp'],'-dbmp','-r300');


t=2;
ind=find(strcmp(risk_type,types{t}));
risks=T1.Risk_factors(ind);pval=T1.p_1tail(ind);
C_perc=T1.Control_perc(ind);F_perc=T1.Frail_perc(ind);
CF_OR=T1.C_F_OR(ind);FC_OR=T1.F_C_OR(ind);
chi=T1.Chi2(ind);
sum_perc=C_perc+F_perc;
for i=1:length(risks)
    if C_perc(i)<F_perc(i)
        sum_perc(i)=-sum_perc(i);
        chi(i)=-chi(i);
    end
end
ind=find(sum_perc<0); C_perc(ind)=-C_perc(ind);F_perc(ind)=-F_perc(ind);
T2=table(risks,C_perc,F_perc,CF_OR,FC_OR,sum_perc,chi,pval);
T2=sortrows(T2,7,'ascend'); x=1:length(risks);y=[T2.C_perc T2.F_perc];
c=customcolormap([0 0.5 1],[0 0 1; 0.8 0.8 0.8; 1 0 0]); k1=length(find(chi<0));
st=floor(128/k1); cmap=c(1:st:129,:);
k2=length(find(chi>0)); st=floor(128/k2);
cmap2=c(129:st:end,:);cmap=[cmap(1:k1,:); cmap2(1:k2,:)];
figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.4, 0.6])
h=superbar(T2.chi,'Orientation','v','BarFaceColor',cmap); box off 
endpt=T2.chi; i=find(endpt>0);endpt(i)=endpt(i)+6;
j=find(endpt<0);endpt(j)=endpt(j)-6; 
OR=[T2.FC_OR(T2.FC_OR>T2.CF_OR); T2.CF_OR(T2.CF_OR>T2.FC_OR)];
OR=round(OR*100)/100; OR(OR==Inf)=1.00;OR=abs(OR);
for i=1:length(T2.risks)
    if T2.pval(i)<0.05 && T2.pval(i)>0.01; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')']; 
        elseif T2.pval(i)<0.01 && T2.pval(i)>0.001; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')'];
        elseif T2.pval(i)<0.001; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')']; 
    else pre=''; post=''; end
    
    if endpt(i)<0; drc='right';drc2='left'; bgn=0.5; else; drc='left'; drc2='right';bgn=-0.5; end
    text(i,bgn,[pre T2.risks{i} ' (' num2str(OR(i)) ')'],'HorizontalAlignment',drc2,'Rotation',40,'fontsize',9);
end
k=round(max(abs(endpt))); ylim([-k k]); %legend({'Control','Cognitively Frail'},'NumColumns',2,'Location','SouthOutside')
grid on; %ylabel('Nutrition, smoking alcohol & use','Fontweight','bold'); 
set(gca,'YColor','none'); set(gca,'XTick',[],'XColor','w')
%xtickangle(90)
print(gcf,['/imaging/ek01/camcan3000_EK/physical_health_bar_plots.bmp'],'-dbmp','-r300');


t=3;
ind=find(strcmp(risk_type,types{t}));
risks=T1.Risk_factors(ind);pval=T1.p_1tail(ind);
C_perc=T1.Control_perc(ind);F_perc=T1.Frail_perc(ind);
CF_OR=T1.C_F_OR(ind);FC_OR=T1.F_C_OR(ind);
chi=T1.Chi2(ind);
sum_perc=C_perc+F_perc;
for i=1:length(risks)
    if C_perc(i)<F_perc(i)
        sum_perc(i)=-sum_perc(i);
        chi(i)=-chi(i);
    end
end
ind=find(sum_perc<0); C_perc(ind)=-C_perc(ind);F_perc(ind)=-F_perc(ind);
T2=table(risks,C_perc,F_perc,CF_OR,FC_OR,sum_perc,chi,pval);
T2=sortrows(T2,7,'ascend'); x=1:length(risks);y=[T2.C_perc T2.F_perc];
c=customcolormap([0 0.5 1],[0 0 1; 0.8 0.8 0.8; 1 0 0]); k1=length(find(chi<0));
st=floor(128/k1); cmap=c(1:st:129,:);
k2=length(find(chi>0)); st=floor(128/k2);
cmap2=c(129:st:end,:);cmap=[cmap(1:k1,:); cmap2(1:k2,:)];

figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.3, 0.5])
h=superbar(T2.chi,'Orientation','v','BarFaceColor',cmap); box off 
endpt=T2.chi; i=find(endpt>0);endpt(i)=endpt(i)+6;
j=find(endpt<0);endpt(j)=endpt(j)-6; 
OR=[T2.FC_OR(T2.FC_OR>T2.CF_OR); T2.CF_OR(T2.CF_OR>T2.FC_OR)];
OR=round(OR*100)/100; OR(OR==Inf)=1.00;OR=abs(OR);
for i=1:length(T2.risks)
    if T2.pval(i)<0.05 && T2.pval(i)>0.01; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')']; 
        elseif T2.pval(i)<0.01 && T2.pval(i)>0.001; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')'];
        elseif T2.pval(i)<0.001; pre='\bf'; post=[' (' num2str(T2.pval(i)) ')']; 
    else pre=''; post=''; end
    
    if endpt(i)<0; drc='right';drc2='left'; bgn=1; else; drc='left'; drc2='right';bgn=-1; end
    text(i,bgn,[pre T2.risks{i} ' (' num2str(OR(i)) ')'],'HorizontalAlignment',drc2,'Rotation',30,'fontsize',9);
end
k=round(max(abs(endpt))); ylim([-k k]); %legend({'Control','Cognitively Frail'},'NumColumns',2,'Location','SouthOutside')
grid on; %ylabel('Nutrition, smoking alcohol & use','Fontweight','bold'); 
set(gca,'YColor','none'); set(gca,'XTick',[],'XColor','w')

print(gcf,['/imaging/ek01/camcan3000_EK/social_demo_bar_plots.bmp'],'-dbmp','-r300');

%% Test continuous measures

mat=readmatrix('cont_var.txt'); c_ind=1:190; f_ind=191:380;
varnm2={'age','yoe','MMSE','ACER','size of household','n social groups','age started smoking','n cigarettes pd','n cooked veg pd','n raw veg pd','n fruits pd','n dry fruits pd','n bread slices pw','n cups of tea pd','n cups of coffee pd','n water pw','n red wine pw','n white wine pw','n beer pw','n spirits pw','n fort wine pw','R ear 1000Hz tones','R ear 3000Hz tones','L ear 1000Hz tones','L ear 3000Hz tones','hypertension age','high cholesterol age','angina age','heart attack age','arryhtmia age','varicose veins age','migraine age','stroke age','pulmonary embolism age','DVT age','vascular disease age','diabetes age','thyroid disease age','peptic ulcer age','polyps age','gallstones age','pancreatitis age','appendicitis age','liver disease age','hayfever eczema age','asthma age','bronchitis emphysema age','allergies age','arthritis age','osteoporosis age','tuberculosis age','enlarged prostate age','insomnia age','depression age','other psychiatric illness age','benign growths age','cancer age','intermittent claudication age','chronic bronchitis age','epilepsy age','anaesthetic age','n anaesthetic','head injury age','hip fracture age','wrist fracture age','vertebrae fracture age','period start age','period stop age','HADS anxiety','HADS depression','retirement age','partner retirement age','n children','PAEE','home PAEE','work PAEE','leisure PAEE','left eye ordinal','right eye ordinal','both eyes ordinal','weekday sleep','weekend sleep','bed time','min to fall asleep','wake up time','hours of sleep'};
% 
% for v=5:length(varnm2)
%     a=mat(:,v);a(find(a==66))=NaN; a(find(a==77))=NaN; a(find(a==88))=NaN; a(find(a==99))=NaN;
%     mat(:,v)=a;
% end; mat(find(mat>100))=NaN;

for v=1:size(mat,2)
    [h p ci stats]=ttest2(mat(c_ind,v),mat(f_ind,v),'vartype','unequal');
    tval(v,1)=stats.tstat; degf(v,1)=stats.df; pval2(v,1)=p;pval1(v,1)=p/2;
    c_m(v,1)=nanmean(mat(c_ind,v)); f_m(v,1)=nanmean(mat(f_ind,v));
    c_sd(v,1)=nanstd(mat(c_ind,v)); f_sd(v,1)=nanstd(mat(f_ind,v));
end

T2=table(varnm2',c_m,c_sd,f_m,f_sd,tval,degf,pval2,pval1,'VariableNames',{'Risk_factors','C_mean','C_SD','F_mean','F_SD','t','df','p_2tail','p_1tail'})
save([outdir 'continuous_risk_factors_ttest_results.mat'],'T2');

T2(find(strcmp(T2.Risk_factors,'MMSE')),:)=[];
T2(find(strcmp(T2.Risk_factors,'ACER')),:)=[];
T2(find(isnan(T2.t)),:)=[];

T2=sortrows(T2,6,'ascend'); 
c=customcolormap([0 0.5 1],[0 0 1; 0.8 0.8 0.8; 1 0 0]); k1=length(find(T2.t<0));
st=floor(128/k1); cmap=c(1:st:129,:);
k2=length(find(T2.t>0))+length(find(T2.t==0)); st=floor(128/k2);
cmap2=c(129:st:end,:);cmap=[cmap(1:k1,:); cmap2(1:k2,:)];

figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.45, 0.5])
h=superbar(T2.t,'Orientation','v','BarFaceColor',cmap); box off 
endpt=T2.t; i=find(endpt>0);endpt(i)=endpt(i)+6;
j=find(endpt<0);endpt(j)=endpt(j)-6; 
for i=1:length(T2.Risk_factors)
    if T2.p_1tail(i)<0.05;  pre='\bf'; post=[' (' num2str(T2.p_1tail(i)) ')']; 
    else pre=''; post=''; end
    
    if endpt(i)<0; drc='right';drc2='left'; bgn=0.2; else; drc='left'; drc2='right';bgn=-0.2; end
    text(i,bgn,[pre T2.Risk_factors{i} ' ' post],'HorizontalAlignment',drc2,'Rotation',30,'fontsize',9);
end
k=round(max(abs(endpt))); ylim([-k k]); %legend({'Control','Cognitively Frail'},'NumColumns',2,'Location','SouthOutside')
grid on; %ylabel('Nutrition, smoking alcohol & use','Fontweight','bold'); 
set(gca,'YColor','none'); set(gca,'XTick',[],'XColor','w')
print(gcf,['/imaging/ek01/camcan3000_EK/cont_variables_bar_plots.bmp'],'-dbmp','-r300');

k2=find(T2.t<0);
figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.25, 0.4])
h=superbar(T2.t(k2),'Orientation','v','BarFaceColor',cmap(k2,:)); box off 
endpt=T2.t(k2); i=find(endpt>0);endpt(i)=endpt(i)+6;
j=find(endpt<0);endpt(j)=endpt(j)-6; 
for i=1:max(k2)
    if T2.p_1tail(i)<0.05;  pre='\bf'; post=[' (' num2str(T2.p_1tail(i)) ')']; 
    else pre=''; post=''; end
    
    if endpt(i)<0; drc='right';drc2='left'; bgn=0.2; else; drc='left'; drc2='right';bgn=-0.2; end
    text(i,bgn,[pre T2.Risk_factors{i} ' ' post],'HorizontalAlignment',drc2,'Rotation',30,'fontsize',9);
end
k=round(max(abs(endpt))); ylim([-4 3]); %legend({'Control','Cognitively Frail'},'NumColumns',2,'Location','SouthOutside')
grid on; %ylabel('Nutrition, smoking alcohol & use','Fontweight','bold'); 
set(gca,'YColor','none'); set(gca,'XTick',[],'XColor','w')
print(gcf,['/imaging/ek01/camcan3000_EK/cont_variables_bar_plots_negative_end.bmp'],'-dbmp','-r300');

k2=[find(T2.t==0); find(T2.t>0)];
figure('color','w'); set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.25, 0.4])
h=superbar(T2.t(k2),'Orientation','v','BarFaceColor',cmap(k2,:)); box off 
endpt=T2.t(k2); i=find(endpt>0);endpt(i)=endpt(i)+6;
j=find(endpt<0);endpt(j)=endpt(j)-6; 
k=1;
for i=min(k2):max(k2)
    if T2.p_1tail(i)<0.05;  pre='\bf'; post=[' (' num2str(T2.p_1tail(i)) ')']; 
    else pre=''; post=''; end
    
    if endpt(k)<0; drc='right';drc2='left'; bgn=0.2; else; drc='left'; drc2='right';bgn=-0.2; end
    text(k,bgn,[pre T2.Risk_factors{i} ' ' post],'HorizontalAlignment',drc2,'Rotation',30,'fontsize',9);
    k=k+1;
end
ylim([-2 5]); %legend({'Control','Cognitively Frail'},'NumColumns',2,'Location','SouthOutside')
grid on; %ylabel('Nutrition, smoking alcohol & use','Fontweight','bold'); 
set(gca,'YColor','none'); set(gca,'XTick',[],'XColor','w')
print(gcf,['/imaging/ek01/camcan3000_EK/cont_variables_bar_plots_positive_end.bmp'],'-dbmp','-r300');

%% Stepwise logistic regression

mat=readmatrix('variables_for_logistic_reg.txt'); c_ind=1:190; f_ind=191:380;

varnm3={'age','income','yoe','n social groups','adult learning','texts relatives','texts friends','phones friends','phones relatives','relatives nearby','has children','rents the house','drinking','white wine','red wine','beer','current smoker','raw veg','poultry','cheese','oily fish','beef','tea','coffee','both eyes ordinal','hearing 1000','mobility issues','SMI','varicose veins','angina','chronic bronchitis','peptic ulcer','int claudication','HADS depression','dental problems','cancer','appendicitis','min to fall asleep','takes sleeping drugs','hrs slept'};
varnm3{end+1}='group';

group=[]; group(c_ind,1)=1; group(f_ind,1)=0;
mat(:,end+1)=group;

for i=1:size(mat,2)
    for j=1:size(mat,2)
        [r p]=corrcoef(mat(:,i),mat(:,j),'rows','complete');
        rmat(i,j)=r(1,2); pmat(i,j)=p(1,2);
    end
end; 

figure('color','w'); a=rmat;imagesc(a);axis equal;colormap('redblue');
caxis([-0.5 0.5]); xlim([0.5 size(mat,2)+.5]); xlim([0.5 size(mat,2)+.5]);
set(gca,'XTick',1:size(mat,2),'YTick',1:size(mat,2),'XTickLabels',varnm3,'YTickLabels',varnm3)
xtickangle(40); colorbar
b=pmat;
[x,y]=ind2sub(size(b),find(b<=0.05)); hold on;
for d=1:length(x)
    plot(y(d),x(d),'-s','markersize',3,'markeredgecolor','k','markerfacecolor','k')
end
print(gcf,['/imaging/ek01/camcan3000_EK/correlation_plot.bmp'],'-dbmp','-r300');

% tree=linkage(rmat,'average');
% figure('color','w');set(gcf, 'Units', 'normal', 'Position', [0, 0, 0.3, 0.4])
% H = dendrogram(tree,0,'ColorThreshold',1.6,'Labels',varnm3);
% set(H,'LineWidth',1.5); xtickangle(40); 

H=clustergram(rmat,'Colormap','redblue','ColumnLabels',varnm3,...
    'RowLabels',varnm3,'ColumnLabelsRotate',40,...
    'DisplayRange',0.3,'Dendrogram',1.6, 'Annotate',false,'DisplayRatio',1/8); %click print to figure
print(gcf,['/imaging/ek01/camcan3000_EK/clustergram.bmp'],'-dbmp','-r300');

vif=diag(inv(rmat));

ind=[2:6 8:40];
mdl=stepwiseglm(mat(:,ind),group,'Constant','upper','linear','VarNames',varnm3([ind 41]))
mdl=stepwiseglm(mat(:,ind),group,'linear','Criterion','rsquared','VarNames',varnm3([ind 41]))
% mdl=stepwiseglm(mat(:,ind),group,'linear','Criterion','bic','VarNames',varnm3([ind 41]))


mat=readmatrix('variables_for_logistic_reg.txt'); 
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan')); % zscoring function that omits nans

for i=1:size(mat,2)
    mat2(:,i)=zscor_xnan(mat(:,i));
end

[coeff score latent t2 explained]=pca(mat2(:,2:40),'rows','complete','NumComponents',15,'Centered',false,'Algorithm','als');
figure; bar(explained)
figure;imagesc(coeff); axis equal; colormap('redblue'); caxis([-0.4 0.4])

mdl=stepwiseglm(score(:,1:15),group,'constant','upper','linear',...
    'VarNames',{'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','group'})

mdl=stepwiseglm(score(:,1:15),group,'linear','Criterion','rsquared',...
    'VarNames',{'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','group'})



