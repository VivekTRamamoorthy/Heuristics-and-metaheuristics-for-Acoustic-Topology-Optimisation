%% GA from benchmark

runfilename=[runfname 'GA'];
% cd GA
% cd(runfilename)
% load(runfilename);
% cd ..
% cd ..
add2Path

%% Memoization
% objfunc=@griddecoder_comb;
objfunc=memoize(@griddecoder_comb);
objfunc.CacheSize=100;
stat=stats(objfunc);

%% Initiating savefile
cd GA
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

%% Genetic Algorithm controls
solutionspace=[0 1];
genmax=(budget-npop)/2;gencounter=0;
crossovertype='uniform'; % 'conservativerandombit';%'conservativesinglepoint';
mutationtype='randomize';%'swap';
mutationrate=1/length(mesh.N);
selectiontype='tournament'; toursize=2;
replacementmethod='weakelitism';

%% Random number seed
rng(123)

res.GA.trial(n_trials).bestgenomes=zeros(genomelength,genmax+1);
res.GA.trial(n_trials).bestfitnesses=zeros(1,genmax+1);
% res.GA.trial(n_trials).popevol=zeros(genomelength,npop,genmax+1);
res.GA.trial(n_trials).genomes=zeros(genomelength,budget);
res.GA.trial(n_trials).fitnesses=zeros(1,budget);
res.GA.trial(n_trials).fitevol=zeros(genmax,npop);
res.GA.avgfitevoltrial=zeros(genmax+1,n_trials);
res.GA.bestfitevoltrial=zeros(genmax+1,n_trials);

GUgenomes=zeros(genomelength,n_trials*budget);GUGcolfilled=0;
GUfitnesses=zeros(1,n_trials*budget);

tic
begintime=toc;

for trial=1:n_trials
    
    %% initializing population
    
    popgenome=zeros(mesh.N,npop);
    popfitness=zeros(1,npop);
    for i=1:npop

        popgenome(:,i)=randi(2,mesh.N,1)-1; % for no vol constraint
        popfitness(1,i)=objfunc(popgenome(:,i),parameters);
        disp(['Initializing ' num2str(i) ' out of ' num2str(npop)])

    end
    
    %% Genetic algorithm
    
    GAcard=['GA:  npop ' num2str(npop) ' gens ' num2str(genmax) ' sel ' selectiontype  '2' ' C '  crossovertype ' M ' mutationtype ' MR ' num2str(mutationrate) ' repl ' replacementmethod ];
    [popgenome,popfitness,bestgenome,bestfitness,fitevol,gencounter,avgfitevol,popevol]=... 
        genetic_algorithm(objfunc,popgenome,popfitness,solutionspace,genmax,selectiontype,crossovertype,mutationtype,mutationrate,replacementmethod,[],gencounter);
    
    %% saving data in structure
    
    res.GA.trial(trial).bestgenomes=bestgenome;
    res.GA.trial(trial).bestfitnesses=bestfitness;
    [res.GA.bestfitnesstrial(1,trial),maxindex]=max(bestfitness);
    res.GA.bestgenometrial(:,trial)=bestgenome(:,maxindex);
    
%     res.GA.trial(trial).popevol=popevol;
    res.GA.trial(trial).fitevol=fitevol;
    res.GA.avgfitevoltrial(:,trial)=mean(fitevol,2); % each column is avg fitness across generation in each trial
    res.GA.avgfitevol=mean(res.GA.avgfitevoltrial(:,1:trial),2); % avg across trials of avg fitness evolution
    
    res.GA.bestfitevoltrial(:,trial)=max(res.GA.trial(trial).fitevol,[],2);
    res.GA.bestfitevol=mean(res.GA.bestfitevoltrial(:,1:trial),2);
    
    [res.GA.trial(trial).bestfitness,maxindex]=max(bestfitness);
    res.GA.trial(trial).bestgenome=bestgenome(:,maxindex);
    
    %% finding unique genomes
    popevolcolumns=zeros(size(popevol,1),size(popevol,2)*size(popevol,3));
    fitevolcolumns=zeros(1,size(popevol,2)*size(popevol,3));
    for i=1:size(popevol,3) %  converting popevol to columns
        columns=(i-1)*size(popevol,2)+1:i*size(popevol,2);
        popevolcolumns(:,columns)=popevol(:,:,i);
        fitevolcolumns(1,columns)=fitevol(i,:);
    end
    [ugenome,uniquecols]=unique(popevolcolumns','rows');
    ugenome=ugenome';
    uniquecols=uniquecols';
    ufitness=fitevolcolumns(uniquecols);
    
    GUgenomes(:,GUGcolfilled+1:GUGcolfilled+length(ufitness))=ugenome;
    GUfitnesses(1,GUGcolfilled+1:GUGcolfilled+length(ufitness))=ufitness;
    GUGcolfilled=GUGcolfilled+length(ufitness);
%     res.GA.trial(trial).ugenome=ugenome;
%     res.GA.trial(trial).ufitness=ufitness;
    
    %% saving mat file after finishing the run
    cd GA
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
    ETC=currenttime/curriteration*(totiteration);
    ETR=ETC-currenttime;
    ETRday=floor(ETR/(24*3600));
    ETRhr=floor(mod(ETR,24*3600)/3600);
    ETRmin=floor(mod(ETR,3600)/60);
    ETRsec=mod(ETR,60);
    disp('--------- Remaining time estimate----------')
    disp(['Avg time per trial: ' num2str(currenttime/curriteration)])
    disp(['Est. time rem: ' num2str(ETRday) ' D ' num2str(ETRhr) ' H ' num2str(ETRmin) ' M ' num2str(ETRsec) ' S'])
    disp(['Percent complete : ' num2str(percent_complete) ' % '])
    disp('-------------------------------------------')
    
end


%% post proc

[res.GA.bestfitness,besttrial]=max(res.GA.bestfitnesstrial);
res.GA.bestgenome=res.GA.bestgenometrial(:,besttrial);

GUgenomes(:,GUGcolfilled+1:end)=[];
GUfitnesses(:,GUGcolfilled+1:end)=[];
    
res.GA.fevalsgen=npop:2:budget;

disp('Hurray! Genetic algorithm successfully completed.')

cd GA
cd(runfilename)
save(runfilename)
cd ..
cd ..

endtime=toc;
disp(['Total time elapsed: ' num2str(endtime)])

%% Plotting average fitness evolution

figure(17)
box on

plot(res.GA.fevalsgen,res.GA.avgfitevol,'b-','linewidth',2,'displayname','Average fitness in population')

hold on
plot(res.GA.fevalsgen,res.GA.bestfitevol,'r--','linewidth',2,'displayname','Best fitness in population')

xlabel('Function evaluations','interpreter','latex')
ylabel('Average fitness across trials','interpreter','latex')
title('GA','interpreter','latex')
ylim([0 1])
legend('location','southeast','interpreter','latex')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])

cd GA
cd(runfilename)
saveas(gcf,[runfilename '_fitevol.fig'],'fig')
saveas(gcf,[runfilename '_fitevol.eps'],'epsc')
saveas(gcf,[runfilename '_fitevol.png'],'png')
cd ..
cd ..

%% plotting fitness evolution for all trials
figure(18)

plot(res.GA.avgfitevoltrial,'b','linewidth',1,'displayname','Average fitness in population')

hold on

plot(res.GA.bestfitevoltrial,'r','linewidth',1,'displayname','Best fitness in population')

xlabel('Generation number','interpreter','latex')
ylabel('Trialwise fitnesses','interpreter','latex')
title('GA: blue-avg red-best','interpreter','latex')
ylim([0 1])

% legend({'Average fitnesses','Best fitnesses'},'interpreter','latex','location','southeast')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])

cd GA
cd(runfilename)
saveas(gcf,[runfilename '_fitevoltrial.fig'],'fig')
saveas(gcf,[runfilename '_fitevoltrial.eps'],'epsc')
saveas(gcf,[runfilename '_fitevoltrial.png'],'png')
cd ..
cd ..


%% best shape
parameters{4}=1;
alpha=griddecoder_cont(res.GA.bestgenome);
parameters{4}=0;
title(['GA best shape: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

cd GA
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

cd GA
cd(runfilename)
saveas(gcf,[runfilename '_fullshape.fig'],'fig')
saveas(gcf,[runfilename '_fullshape.eps'],'epsc')
saveas(gcf,[runfilename '_fullshape.png'],'png')
cd ..
cd ..
