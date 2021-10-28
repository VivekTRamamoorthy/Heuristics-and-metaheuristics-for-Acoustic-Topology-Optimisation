%% DE from benchmark

runfilename=[runfname 'DE'];
% load(runfilename);
add2Path
mesh.penalty=3;
%% Memoization
objfunc=@griddecoder_cont;
% objfunc=memoize(@griddecoder_cont);
% objfunc.CacheSize=100;
% stat=stats(objfunc);

%% Initiating savefile
cd DE
mkdir(runfilename)
cd(runfilename)
if isfile([runfilename '.mat'])
    warning('A save file with the same name already exists. Rename or create backup') 
    keyboard
else
save(runfilename)
end
cd ..
cd ..

%% Differential evolution controls
solutionspace=[0 1];
genmax=round((budget-npop)/npop);
gencounter=0;
F=0.2;
CR=0.2;
range=[0 1];

%% Random number seed
rng(123)

res.DE.trial(n_trials).bestgenomes=zeros(genomelength,genmax+1);
res.DE.trial(n_trials).bestfitnesses=zeros(1,genmax+1);
% res.DE.trial(n_trials).popevol=zeros(genomelength,npop,genmax+1);
res.DE.trial(n_trials).fitevol=zeros(genmax,npop);
res.DE.avgfitevoltrial=zeros(genmax+1,n_trials);
res.DE.bestfitevoltrial=zeros(genmax+1,n_trials);

tic
begintime=toc;

for trial=1:n_trials
    
    %% initializing population
    
    popgenome=zeros(mesh.N,npop);
    popfitness=zeros(1,npop);
    for i=1:npop
        popgenome(:,i)=rand(mesh.N,1);
        popfitness(1,i)=objfunc(popgenome(:,i),parameters);
        disp(['Initializing ' num2str(i) ' out of ' num2str(npop)])

    end
    
    %% Differential Evolution
    disp('Differential Evolution started')

    [popgenome,popfitness,bestgenome,bestfitness,fitevol,gencounter,avgfitevol,popevol]=...
        differential_evolution(objfunc,popgenome,popfitness,range,genmax,F,CR,gencounter);
    
    %% saving data in structure

    res.DE.trial(trial).bestgenomes=bestgenome; % bestgenome in generation
    res.DE.trial(trial).bestfitnesses=bestfitness; % best fitnesses in generation
    [res.DE.bestfitnesstrial(1,trial),maxindex]=max(bestfitness);
    res.DE.bestgenometrial(:,trial)=bestgenome(:,maxindex);
    
%     res.DE.trial(trial).popevol=popevol;
    res.DE.trial(trial).fitevol=fitevol;
    res.DE.avgfitevoltrial(:,trial)=mean(fitevol,2); % each column is avg fitness across generation in each trial
    res.DE.avgfitevol=mean(res.DE.avgfitevoltrial(:,1:trial),2); % avg across trials of avg fitness evolution
    
    res.DE.bestfitevoltrial(:,trial)=max(res.DE.trial(trial).fitevol,[],2);
    res.DE.bestfitevol=mean(res.DE.bestfitevoltrial(:,1:trial),2);
    
    [res.DE.trial(trial).bestfitness,maxindex]=max(bestfitness);
    res.DE.trial(trial).bestgenome=bestgenome(:,maxindex);
     
    
    
    
    
    %% saving mat file after finishing the run
    cd DE
    cd(runfilename)
    save(runfilename)
    cd ..
    cd ..
    disp(['Trial ' num2str(trial) ' finished'])
    
    
    
    %% calculating percent complete and time remaining
    curriteration=trial;
    totiteration=n_trials;
    percent_complete=curriteration/(totiteration)*100;
    
    currenttime=toc;
    ETR=currenttime/curriteration*(totiteration)-currenttime;
    ETRday=floor(ETR/(24*3600));
    ETRhr=floor(mod(ETR,24*3600)/3600);
    ETRmin=floor(mod(ETR,3600)/60);
    ETRsec=mod(ETR,60);
    disp(['------------- Remaining time estimate----------'])
    disp(['Avg time per trial: ' num2str(currenttime/curriteration)])
    disp(['Est. time remaining: ' num2str(ETRday) ' D ' num2str(ETRhr) ' H ' num2str(ETRmin) ' M ' num2str(ETRsec) ' S'])
    disp(['Percent complete : ' num2str(percent_complete) ' % '])
    disp(['-----------------------------------------------'])
end


%% post proc

[res.DE.bestfitness,besttrial]=max(res.DE.bestfitnesstrial);
res.DE.bestgenome=res.DE.bestgenometrial(:,besttrial);

res.DE.fevalsgen=[npop:npop:(genmax+1)*npop]';

disp('Hurray! Differential evolution successfully completed.')

%% generating the discretized best genomes and shapes

res.DE.bestgenomedtrial=round(res.DE.bestgenometrial);
for trial=1:n_trials
res.DE.bestfitnessdtrial(1,trial)=griddecoder_comb(res.DE.bestgenomedtrial(:,trial));
disp(['trial ' num2str(trial) ' discretized! '])
end

[res.DE.bestfitnessd,besttriald]=max(res.DE.bestfitnessdtrial);
res.DE.bestgenomed=res.DE.bestgenomedtrial(:,besttriald);



cd DE
cd(runfilename)
save(runfilename)
cd ..
cd ..

endtime=toc;
disp(['Total time elapsed: ' num2str(endtime)])

%% plotting threshold cutoff solutions and their fitnesses
% % meshplotter(parameters{1},bestgenome(:,end))
% % temp=popgenome(:,1);
% for cutoff=0:10:100
% % temp=bestgenome(:,end);
% [~,maxindex]=max(popfitness);
% temp=popgenome(:,maxindex);
% [sorted,index]=sort(temp);
% % temp=popgenome(:,10);
% % cutoff=0;
% temp(index(1:cutoff))=0;
% temp(index(cutoff+1:100))=1;
% % meshplotter(parameters{1},popgenome(:,10))
% % meshplotter(parameters{1},temp)
%
% parameters{4}=1;
% objfunc(temp)
% cd DE
% saveas(gcf,['DE_cutoff_' num2str(cutoff) '.fig'],'fig')
% saveas(gcf,['DE_cutoff_' num2str(cutoff) '.eps'],'epsc')
% saveas(gcf,['DE_cutoff_' num2str(cutoff) '.png'],'png')
% cd ..
% end

%% Plotting average fitness evolution

figure(17)
box on

plot(res.DE.fevalsgen,res.DE.avgfitevol,'b-','linewidth',2,'displayname','Average fitness in population')

hold on
plot(res.DE.fevalsgen,res.DE.bestfitevol,'r--','linewidth',2,'displayname','Best fitness in population')

xlabel('Function evaluations','interpreter','latex')
ylabel('Average fitness across trials','interpreter','latex')
title('DE','interpreter','latex')
ylim([0 1])
legend('location','southeast','interpreter','latex')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf ,'position',[440   378   560   420],'Color',[1 1 1])

cd DE
cd(runfilename)
saveas(gcf,[runfilename '_fitevol.fig'],'fig')
saveas(gcf,[runfilename '_fitevol.eps'],'epsc')
saveas(gcf,[runfilename '_fitevol.png'],'png')
cd ..
cd ..

%% plotting fitness evolution for all trials
figure(18)

plot(res.DE.avgfitevoltrial,'b','linewidth',1,'displayname','Average fitness in population')

hold on

plot(res.DE.bestfitevoltrial,'r','linewidth',1,'displayname','Best fitness in population')

xlabel('Generation number','interpreter','latex')
ylabel('Trialwise fitnesses','interpreter','latex')
title('DE: blue-avg red-best','interpreter','latex')
ylim([0 1])

% legend({'Average fitnesses','Best fitnesses'},'interpreter','latex','location','southeast')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])

cd DE
cd(runfilename)
saveas(gcf,[runfilename '_fitevoltrial.fig'],'fig')
saveas(gcf,[runfilename '_fitevoltrial.eps'],'epsc')
saveas(gcf,[runfilename '_fitevoltrial.png'],'png')
cd ..
cd ..


%% best shape
parameters{4}=1;
alpha=griddecoder_cont(res.DE.bestgenome);
parameters{4}=0;
title(['DE best shape: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

cd DE
cd(runfilename)
saveas(gcf,[runfilename '_bestshape.fig'],'fig')
saveas(gcf,[runfilename '_bestshape.eps'],'epsc')
saveas(gcf,[runfilename '_bestshape.png'],'png')
cd ..
cd ..

%% flatshape
parameters{4}=1;
alpha=griddecoder_cont(ones(size(bestgenome(:,end))));
parameters{4}=0;
title(['Flat layer: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

cd DE
cd(runfilename)
saveas(gcf,[runfilename '_fullshape.fig'],'fig')
saveas(gcf,[runfilename '_fullshape.eps'],'epsc')
saveas(gcf,[runfilename '_fullshape.png'],'png')
cd ..
cd ..
