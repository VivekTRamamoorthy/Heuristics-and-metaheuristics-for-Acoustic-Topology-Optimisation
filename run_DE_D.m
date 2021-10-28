%% DE_D from benchmark
algorithm='DE_D';
runfilename=[runfname 'DE_D'];
% load(runfilename);
add2Path
mesh.penalty=3;
%% Memoization
% objfunc=@griddecoder_cont;
objfunc=memoize(@griddecoder_comb);
objfunc.CacheSize=100;
stat=stats(objfunc);

%% Initiating savefile 
cd DE_D
mkdir(runfilename)
cd(runfilename)
if isfile([runfilename '.mat'])
    warning('A save file with the same name already exists. Rename or create backup!')
    %keyboard
    mode= input('Enter, y to continue replacing, r to restart, c to cancel \n','s');
    switch mode
        case 'y'
            trialbegin=1;

            disp('Are you sure you want to continue replacing existing results?')
            mode='fresh';
            disp('change the trialvector here if you want specific trials to run in this machine');
            disp('if you dont change trial vector, all trials will run')
            trialvector=trialbegin:n_trials;
            keyboard
        case 'c'
            disp('Cancelled run')
            cd ..
            cd ..
            return
        case 'r'
            load(runfilename)
            if ~isnumeric('trial')
                trial=0;
            end
            trialbegin=trial+1;
            disp(['Current trial in the existing file is ' num2str(trial)])
            disp('Are you sure you want to restart?')
            disp('Change the trialvector here if you want specific trials to run in this machine.');
            % all remaining trials
            trialvector=trialbegin:n_trials;
            
            keyboard
            % trialvector=6:10;
            
        otherwise
            warning('invalid input')
            disp('Continue replacing the existing results? ')
            disp('specify the trial number to be run on this machine in trialvector');
            disp('if you dont change trial vector, all trials will be run')
            mode='fresh';
            trialbegin=1;
            trialvector=trialbegin:n_trials;
            keyboard
    end
else
    mode='fresh';
    trialbegin=1;
    trialvector=trialbegin:n_trials;
    save(runfilename)
end
cd ..
cd ..

%% setting inputs for fresh runs

if strcmp(mode,'fresh')
    %% Differential evolution controls
    solutionspace=[0 1];
    genmax=round((budget-npop)/npop);
    gencounter=0;
    F=0.8;
    CR=0.9;
    range=[0 1];
    
    %% Initialization new datastructure for results
    
    res.DE_D.trial(n_trials).bestgenomes=zeros(genomelength,genmax+1);
    res.DE_D.trial(n_trials).bestfitnesses=zeros(1,genmax+1);
    % res.DE_D.trial(n_trials).popevol=zeros(genomelength,npop,genmax+1);
    res.DE_D.trial(n_trials).fitevol=zeros(genmax,npop);
    res.DE_D.avgfitevoltrial=zeros(genmax+1,n_trials);
    res.DE_D.bestfitevoltrial=zeros(genmax+1,n_trials);
    
    trialscomplete=zeros(n_trials,1);
end

%% Run from here for restarts
tic
begintime=toc;

for trial=trialvector
    disp(['Trial ' num2str(trial) ' initiated. Sit back and relax.'])
    %% Trial based Random number seed
    disp('Using random number seed based on trial number')
    rngseed=123+(trial-1)*1000000;
    rng(rngseed)
    %% initializing population
    
    popgenome=zeros(mesh.N,npop);
    popfitness=zeros(1,npop);
    for i=1:npop
        popgenome(:,i)=rand(mesh.N,1);
        popfitness(1,i)=objfunc(round(popgenome(:,i)),parameters);
        disp(['Initializing ' num2str(i) ' out of ' num2str(npop)])
        
    end
    
    %% Differential Evolution
    disp('Differential Evolution started')
    
    [popgenome,popfitness,bestgenome,bestfitness,fitevol,gencounter,avgfitevol,popevol]=...
        differential_evolution_discrete(objfunc,popgenome,popfitness,range,genmax,F,CR,gencounter);
    % differential_evolution_discrete
    %% saving data in structure
    
    res.DE_D.trial(trial).bestgenomes=bestgenome; % bestgenome in generation
    res.DE_D.trial(trial).bestfitnesses=bestfitness; % best fitnesses in generation
    [res.DE_D.bestfitnesstrial(1,trial),maxindex]=max(bestfitness);
    res.DE_D.bestgenometrial(:,trial)=bestgenome(:,maxindex);
    
    %     res.DE_D.trial(trial).popevol=popevol;
    res.DE_D.trial(trial).fitevol=fitevol;
    res.DE_D.avgfitevoltrial(:,trial)=mean(fitevol,2); % each column is avg fitness across generation in each trial
    res.DE_D.avgfitevol=mean(res.DE_D.avgfitevoltrial(:,1:trial),2); % avg across trials of avg fitness evolution
    
    res.DE_D.bestfitevoltrial(:,trial)=max(res.DE_D.trial(trial).fitevol,[],2);
    res.DE_D.bestfitevol=mean(res.DE_D.bestfitevoltrial(:,1:trial),2);
    
    [res.DE_D.trial(trial).bestfitness,maxindex]=max(bestfitness);
    res.DE_D.trial(trial).bestgenome=bestgenome(:,maxindex);
    
    
    
    trialscomplete(trial)=1;
    
    %% saving mat file after finishing the run
    cd DE_D
    cd(runfilename)
    save(runfilename)
    %if mode=='r'
        trialdata=res.DE_D.trial(trial);
        save([runfilename '_trial_' num2str(trial)],'trialdata')
    %end
    cd ..
    cd ..
    disp(['Trial ' num2str(trial) ' finished'])
    
    
    
    %% calculating percent complete and time remaining
    curriteration=trial-trialbegin+1;
    totiteration=n_trials-trialbegin+1;
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

[res.DE_D.bestfitness,besttrial]=max(res.DE_D.bestfitnesstrial);
res.DE_D.bestgenome=res.DE_D.bestgenometrial(:,besttrial);

res.DE_D.fevalsgen=[npop:npop:(genmax+1)*npop]';

disp('Hurray! Differential evolution successfully completed.')

%% generating the discretized best genomes and shapes
% 
res.DE_D.bestgenomedtrial=round(res.DE_D.bestgenometrial);

% for trial=1:n_trials
%     res.DE_D.bestfitnessdtrial(1,trial)=griddecoder_comb(res.DE_D.bestgenomedtrial(:,trial));
%     disp(['trial ' num2str(trial) ' discretized! '])
% end

[res.DE_D.bestfitnessd,besttriald]=max(res.DE_D.bestfitnesstrial);
res.DE_D.bestgenomed=res.DE_D.bestgenomedtrial(:,besttriald);

%%

cd DE_D
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
% cd DE_D
% saveas(gcf,['DE_D_cutoff_' num2str(cutoff) '.fig'],'fig')
% saveas(gcf,['DE_D_cutoff_' num2str(cutoff) '.eps'],'epsc')
% saveas(gcf,['DE_D_cutoff_' num2str(cutoff) '.png'],'png')
% cd ..
% end

%% Plotting average fitness evolution

figure(17)
box on

plot(res.DE_D.fevalsgen,res.DE_D.avgfitevol,'b-','linewidth',2,'displayname','Average fitness in population')

hold on
plot(res.DE_D.fevalsgen,res.DE_D.bestfitevol,'r--','linewidth',2,'displayname','Best fitness in population')

xlabel('Function evaluations','interpreter','latex')
ylabel('Average fitness across trials','interpreter','latex')
title('DE_D','interpreter','latex')
ylim([0 1])
legend('location','southeast','interpreter','latex')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf ,'position',[50   50   560   420],'Color',[1 1 1])

cd DE_D
cd(runfilename)
saveas(gcf,[runfilename '_fitevol.fig'],'fig')
saveas(gcf,[runfilename '_fitevol.eps'],'epsc')
saveas(gcf,[runfilename '_fitevol.png'],'png')
cd ..
cd ..

%% plotting fitness evolution for all trials
figure(18)

plot(res.DE_D.avgfitevoltrial,'b','linewidth',1,'displayname','Average fitness in population')

hold on

plot(res.DE_D.bestfitevoltrial,'r','linewidth',1,'displayname','Best fitness in population')

xlabel('Generation number','interpreter','latex')
ylabel('Trialwise fitnesses','interpreter','latex')
title('DE_D: blue-avg red-best','interpreter','latex')
ylim([0 1])

% legend({'Average fitnesses','Best fitnesses'},'interpreter','latex','location','southeast')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])

cd DE_D
cd(runfilename)
saveas(gcf,[runfilename '_fitevoltrial.fig'],'fig')
saveas(gcf,[runfilename '_fitevoltrial.eps'],'epsc')
saveas(gcf,[runfilename '_fitevoltrial.png'],'png')
cd ..
cd ..


%% best shape
parameters{4}=1;
% alpha=griddecoder_comb(res.DE_D.bestgenome);
meshplotter(mesh,res.DE_D.bestgenomed);
alpha=res.DE_D.bestfitness;
parameters{4}=0;
title(['DE_D best shape: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

cd DE_D
cd(runfilename)
saveas(gcf,[runfilename '_bestshape.fig'],'fig')
saveas(gcf,[runfilename '_bestshape.eps'],'epsc')
saveas(gcf,[runfilename '_bestshape.png'],'png')
cd ..
cd ..
disp('best shape is printed')

%% flatshape
% parameters{4}=1;
% alpha=griddecoder_comb(ones(size(bestgenome(:,end))));
% parameters{4}=0;
% title(['Flat layer: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)
% 
% cd DE_D
% cd(runfilename)
% saveas(gcf,[runfilename '_fullshape.fig'],'fig')
% saveas(gcf,[runfilename '_fullshape.eps'],'epsc')
% saveas(gcf,[runfilename '_fullshape.png'],'png')
% cd ..
% cd ..
