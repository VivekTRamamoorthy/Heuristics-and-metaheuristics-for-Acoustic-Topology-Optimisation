%% run_HC
%% Idea 1a 2d grid optimisation symmetric
algorithm='HC';
runfilename=[runfname 'HC'];
add2Path
%% MEMOIZATION
% objfunc=@griddecoder_cont;
objfunc=memoize(@griddecoder_comb);
objfunc.CacheSize=100;
stat=stats(objfunc);

%% Initiating savefile 
mkdir(algorithm)
cd(algorithm)
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
            mode='r';
            keyboard
            
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
    %% hill climbing
    
    genmax=budget;
    move.type='davisbit'; %'randombitflip' or 'swap'
    acceptance.type='greedy'; %'greedy' or 'simulatedannealing'
    solutionspace=[0 1];
    
    
    %% Initialization new datastructure for results
    
    %     res.HC.trial(n_trials).bestgenomes=zeros(genomelength,genmax);
    %     res.HC.trial(n_trials).bestfitnesses=zeros(1,genmax);
    %     res.HC.trial(n_trials).popevol=zeros(genomelength,npop,genmax+1);
    res.HC.trial(n_trials).fitevol=zeros(genmax,npop);
    res.HC.avgfitevoltrial=zeros(genmax,n_trials);
    res.HC.bestfitevoltrial=zeros(genmax,n_trials);
    
    trialscomplete=zeros(n_trials,1);
	res.HC.fevalsgen=[1:genmax]';
end

%% Run from here for restarts
tic
begintime=toc;

for trial=trialvector
    disp(['Trial ' num2str(trial) ' initiated. Sit back and relax.'])
    %% Trial-based Random number seed
    disp('Using random number seed based on trial number')
    rngseed=123+(trial-1)*1000000;
    rng(rngseed)
    
    %% initializing population
    gencounter=0;
    initialsol=zeros(genomelength,1);
    if VF<1
    initialsol(randperm(genomelength,round(VF*genomelength)))=1;
    elseif VF==1
        initialsol=randi(2,genomelength,1)-1;
        if trial==1
            initialsol=ones(genomelength,1);
        end
    else
        initialsol=randi(2,genomelength,1)-1;
    end
    
    %% Hill climbing
    [bestgenome,bestfitness,fitevol,genomeevol]=hillclimbing...
        (objfunc,initialsol,[],move,acceptance,genmax);
    
    %% saving data in structure
    
    %     res.HC.trial(trial).genomes=genomeevol; % bestgenome in generation
    %     res.HC.trial(trial).bestfitnesses=fitevol; % best fitnesses in generation
    [res.HC.bestfitnesstrial(1,trial),maxindex]=max(fitevol);
    res.HC.bestgenometrial(:,trial)=genomeevol(:,maxindex);
    
    %     res.HC.trial(trial).popevol=popevol;
    res.HC.trial(trial).fitevol=fitevol';
    res.HC.avgfitevoltrial(:,trial)=mean(res.HC.trial(trial).fitevol,2); % each column is avg fitness across generation in each trial
    res.HC.avgfitevol=mean(res.HC.avgfitevoltrial(:,1:trial),2); % avg across trials of avg fitness evolution
    
    res.HC.bestfitevoltrial(:,trial)=cumulativemaximum(max(res.HC.trial(trial).fitevol,[],2));% maximum of the population
    res.HC.bestfitevol=mean(res.HC.bestfitevoltrial(:,1:trial),2);
    
    [res.HC.trial(trial).bestfitness,maxindex]=max(fitevol);
    res.HC.trial(trial).bestgenome=genomeevol(:,maxindex);

    trialscomplete(trial)=1;
    
    %% saving mat file after finishing the trial
    cd(algorithm)
    cd(runfilename)
    % saving all the data
    save(runfilename)
    % saving the trial data
    %if mode=='r'
        trialdata=res.HC.trial(trial);
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

[res.HC.bestfitness,besttrial]=max(res.HC.bestfitnesstrial);
res.HC.bestgenome=res.HC.bestgenometrial(:,besttrial);
res.HC.bestfitevol=cumulativemaximum(res.HC.bestfitevol);

% 	res.HC.fevalsgen=[1:genmax]';

disp('Hurray! Hill climbing successfully completed.')


%% saving runfilename
cd(algorithm)
cd(runfilename)
save(runfilename)
cd ..
cd ..

endtime=toc;
disp(['Total time elapsed: ' num2str(endtime)])
disp('Saved data in:')
disp(runfilename)
% return
%% Plotting average fitness evolution

figure(17)
box on

plot(res.HC.fevalsgen,res.HC.avgfitevol,'b-','linewidth',2,'displayname','Average fitness in population')

hold on
res.HC.bestfitevol=cumulativemaximum(res.HC.bestfitevol);
plot(res.HC.fevalsgen,res.HC.bestfitevol,'r--','linewidth',2,'displayname','Best fitness in population')

xlabel('Function evaluations','interpreter','latex')
ylabel('Average fitness across trials','interpreter','latex')
title('HC','interpreter','latex')
ylim([0 1])
legend('location','southeast','interpreter','latex')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf ,'position',[50   50   560   420],'Color',[1 1 1])

cd HC
cd(runfilename)
saveas(gcf,[runfilename '_fitevol.fig'],'fig')
saveas(gcf,[runfilename '_fitevol.eps'],'epsc')
saveas(gcf,[runfilename '_fitevol.png'],'png')
cd ..
cd ..


%% plotting fitness evolution for all trials
figure(18)

% plot(res.HC.avgfitevoltrial,'b','linewidth',1,'displayname','Average fitness in population')

hold on

plot(res.HC.bestfitevoltrial,'r','linewidth',1,'displayname','Best fitness in population')

xlabel('Generation number','interpreter','latex')
ylabel('Trialwise fitnesses','interpreter','latex')
title('HC: blue-avg red-best','interpreter','latex')
ylim([0 1])

% legend({'Average fitnesses','Best fitnesses'},'interpreter','latex','location','southeast')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])
box on 

cd HC
cd(runfilename)
saveas(gcf,[runfilename '_fitevoltrial.fig'],'fig')
saveas(gcf,[runfilename '_fitevoltrial.eps'],'epsc')
saveas(gcf,[runfilename '_fitevoltrial.png'],'png')
cd ..
cd ..

%% best shape
parameters{4}=1;
% alpha=griddecoder_comb(res.HC.bestgenome);
meshplotter(mesh,res.HC.bestgenome);
alpha=res.HC.bestfitness;
parameters{4}=0;
title(['HC best shape: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)


cd HC
cd(runfilename)
saveas(gcf,[runfilename '_bestshape.fig'],'fig')
saveas(gcf,[runfilename '_bestshape.eps'],'epsc')
saveas(gcf,[runfilename '_bestshape.png'],'png')
cd ..
cd ..
disp('best shape is printed')

%% Completion, recording time and displaying memoize cache


% disp(runfilename)
% stat=stats(objfunc);
% disp(stat)
% stat.Cache.Inputs{1,1}{1,1}
% disp(stat.MostHitCachedInput)
% meshplotter(mesh,stat.MostHitCachedInput.Input{1})









