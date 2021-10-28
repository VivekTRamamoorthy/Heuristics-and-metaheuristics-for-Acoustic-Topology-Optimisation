%% CMA-ES from benchmark

runfilename=[runfname 'CMAES'];
% load(runfilename);
mesh.penalty=3;

%% Memoization
objfunc=@griddecoder_cont; %(maximisation)
% objfunc=memoize(@griddecoder_cont);
% objfunc.CacheSize=100;
% stat=stats(objfunc);

%% Initiating savefile
cd CMAES
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
    
%% CMAES controls
N = mesh.N;           % number of objective variables/problem dimension
npop = 4+floor(3*log(N));  % population size, offspring number
genomelength=N;   
genmax=ceil(budget/npop);



%% Random number seed



for trial=1:n_trials
res.CMAES.trial(trial).bestgenomes=zeros(genomelength,genmax);
res.CMAES.trial(trial).bestfitnesses=zeros(1,genmax);
% res.CMAES.trial(n_trials).popevol=zeros(genomelength,npop,genmax);
res.CMAES.trial(trial).fitevol=zeros(genmax,npop);
end

res.CMAES.avgfitevoltrial=zeros(genmax,n_trials);
res.CMAES.bestfitevoltrial=zeros(genmax,n_trials);

tic
begintime=toc;

for trial=1:n_trials
    %% Trial based Random number seed
    disp('Using random number seed based on trial number')
    rngseed=123+(trial-1)*1000000;
    rng(rngseed)
    
    %% CMAES controls
    
    xmean = rand(N,1);    % objective variables initial point
    sigma = 0.3;          % coordinate wise standard deviation (step size)
    stopcost = 1e-6;   % stop if fitness < stopfitness (minimization)
    stopeval = budget;    % stop after stopeval number of function evaluations

    %% Strategy parameter setting: Selection
    mu = npop/2;               % number of parents/points for recombination
    weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
    mu = floor(mu);
    weights = weights/sum(weights);     % normalize recombination weights array
    mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i
%     genmax=ceil(budget/npop);
    %% Strategy parameter setting: Adaptation
    cc = (4+mueff/N) / (N+4 + 2*mueff/N);  % time constant for cumulation for C
    cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
    c1 = 2 / ((N+1.3)^2+mueff);    % learning rate for rank-one update of C
    cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));  % and for rank-mu update
    damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma
    % usually close to 1
    
    %% Initialize dynamic (internal) strategy parameters and constants
    pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
    B = eye(N,N);                       % B defines the coordinate system
    D = ones(N,1);                      % diagonal D defines the scaling
    C = B * diag(D.^2) * B';            % covariance matrix C
    invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
    eigeneval = 0;                      % track update of B and D
    chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of
    %   ||N(0,I)|| == norm(randn(N,1))
    
    %% Generation Loop
    g=0;
    fevalcounter = 0;  % the next 40 lines contain the 20 lines of interesting code
    
    popgenome=zeros(N,npop);
    popcost=ones(1,npop);
    popfitness=zeros(1,npop);
    
    disp('CMAES started')
    
    while fevalcounter < stopeval
        g=g+1;
        disp(['CMAES generation:  ' num2str(g) ' initiated' ])
        % Generate and evaluate lambda offspring
        for k=1:npop
            popgenome(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C)
            popgenome((popgenome(:,k)>1),k)=1; % constraints
            popgenome((popgenome(:,k)<0),k)=0; % constraints
            
            popfitness(1,k)=objfunc(popgenome(:,k)); % objective function call
            popcost(1,k) = 1-popfitness(1,k);
            
            fevalcounter = fevalcounter+1;
        end

        %% plotting 3D graph
        %       plot3(popgenome(1,:),popgenome(2,:),popfitness(:),'x')
        %       X=linspace(-100,100);
        %       Y=linspace(-50,500);
        %       [X,Y]=meshgrid(X,Y);Z=X;
        %       for i=1:size(X,1)
        %           for j=1:size(X,2)
        %               Z(i,j)=frosenbrock([X(i,j);Y(i,j)]);
        %           end
        %       end
        %       contour(X,Y,Z)
        % Sort by fitness and compute weighted mean into xmean
        [popcost, arindex] = sort(popcost); % minimization sort ascending
        popfitness=popfitness(arindex); % sorting fitness in descending order
        % whereas popgenome is not
        xold = xmean;
        xmean = popgenome(:,arindex(1:mu))*weights;   % recombination, new mean value
        
        % Cumulation: Update evolution paths
        ps = (1-cs)*ps ...
            + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*fevalcounter/npop))/chiN < 1.4 + 2/(N+1);
        pc = (1-cc)*pc ...
            + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
        
        % Adapt covariance matrix C
        artmp = (1/sigma) * (popgenome(:,arindex(1:mu))-repmat(xold,1,mu));
        C = (1-c1-cmu) * C ...                  % regard old matrix
            + c1 * (pc*pc' ...                 % plus rank one update
            + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
            + cmu * artmp * diag(weights) * artmp'; % plus rank mu update
        
        % Adapt step size sigma
        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
        
        % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
        if fevalcounter - eigeneval > npop/(c1+cmu)/N/10  % to achieve O(N^2)
            eigeneval = fevalcounter;
            C = triu(C) + triu(C,1)'; % enforce symmetry
            [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
            D = sqrt(diag(D));        % D is a vector of standard deviations now
            invsqrtC = B * diag(D.^-1) * B';
        end
        
        % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable
        if popcost(1) <= stopcost || max(D) > 1e7 * min(D)
            disp('Fitness goal was reached before end of generations')

            for gi=g:genmax
                res.CMAES.trial(trial).bestgenomes(:,gi) = popgenome(:, arindex(1)); % bestgenome
                res.CMAES.trial(trial).bestfitnesses(1,gi)=popfitness(1,arindex(1)); % bestfitness
                
                res.CMAES.trial(trial).fitevol(gi,:)=popfitness;
            end
            break;
        end
        
        res.CMAES.trial(trial).bestgenomes(:,g) = popgenome(:, arindex(1)); % bestgenome % since popgenome is not sorted in ascending order
        %  res.CMAES.trial(trial).bestfitnesses(1,g)=popfitness(1,arindex(1)); % bestfitness
        res.CMAES.trial(trial).bestfitnesses(1,g)=popfitness(1,1); % bestfitness : popfitness is sorted

        
        %         res.CMAES.trial(trial).popevol(:,:,g)=popgenome;
        res.CMAES.trial(trial).fitevol(g,:)=popfitness;

    
    
    
    
    end % end of generation loop
    
 
    
    %% saving data in structure
    
    %     res.CMAES.trial(trial).bestgenomes=bestgenome;
    %     res.CMAES.trial(trial).bestfitnesses=bestfitness;
    [res.CMAES.bestfitnesstrial(1,trial),maxindex]=max(res.CMAES.trial(trial).bestfitnesses);
    res.CMAES.bestgenometrial(:,trial)=res.CMAES.trial(trial).bestgenomes(:,maxindex);
    
    %     res.CMAES.trial(trial).popevol=popevol;
    %     res.CMAES.trial(trial).fitevol=fitevol;
    res.CMAES.avgfitevoltrial(:,trial)=mean(res.CMAES.trial(trial).fitevol,2); % each column is avg fitness across generation in each trial
    res.CMAES.avgfitevol=mean(res.CMAES.avgfitevoltrial(:,1:trial),2); % avg across trials of avg fitness evolution
    
    res.CMAES.bestfitevoltrial(:,trial)=max(res.CMAES.trial(trial).fitevol,[],2);
    res.CMAES.bestfitevol=mean(res.CMAES.bestfitevoltrial(:,1:trial),2);
    
    [res.CMAES.trial(trial).bestfitness,maxindex]=max(res.CMAES.trial(trial).bestfitnesses);
    res.CMAES.trial(trial).bestgenome=res.CMAES.trial(trial).bestgenomes(:,maxindex);
    
    %% saving mat file after finishing the run
    cd CMAES
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

[res.CMAES.bestfitness,besttrial]=max(res.CMAES.bestfitnesstrial);
res.CMAES.bestgenome=res.CMAES.bestgenometrial(:,besttrial);

res.CMAES.fevalsgen=npop:npop:genmax*npop;
disp('Hurray! CMA Evolution Strategy successfully completed.')

%% discretizing best genomes and shapes

res.CMAES.bestgenomedtrial=round(res.CMAES.bestgenometrial);
for trial=1:n_trials
res.CMAES.bestfitnessdtrial(1,trial)=griddecoder_comb(res.CMAES.bestgenomedtrial(:,trial));
disp(['trial ' num2str(trial) ' discretized! '])
end

[res.CMAES.bestfitnessd,besttriald]=max(res.CMAES.bestfitnessdtrial);
res.CMAES.bestgenomed=res.CMAES.bestgenomedtrial(:,besttriald);


endtime=toc;

cd CMAES
cd(runfilename)
save(runfilename)
cd ..
cd ..


disp(['Total time elapsed: ' num2str(endtime)])

%% Plotting average fitness evolution

figure(17)
box on

plot(npop:npop:genmax*npop,res.CMAES.avgfitevol,'b-','linewidth',2,'displayname','Average fitness in population')

hold on
plot(npop:npop:genmax*npop,res.CMAES.bestfitevol,'r--','linewidth',2,'displayname','Best fitness in population')

xlabel('Function evaluations','interpreter','latex')
ylabel('Average fitness across trials','interpreter','latex')
title('CMAES','interpreter','latex')
ylim([0 1])
legend('location','southeast','interpreter','latex')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])

cd CMAES
cd(runfilename)
saveas(gcf,[runfilename '_fitevol.fig'],'fig')
saveas(gcf,[runfilename '_fitevol.eps'],'epsc')
saveas(gcf,[runfilename '_fitevol.png'],'png')
cd ..
cd ..

%% plotting fitness evolution for all trials
figure(18)

plot(res.CMAES.avgfitevoltrial,'b','linewidth',1,'displayname','Average fitness in population')

hold on

plot(res.CMAES.bestfitevoltrial,'r','linewidth',1,'displayname','Best fitness in population')

xlabel('Generation number','interpreter','latex')
ylabel('Trialwise fitnesses','interpreter','latex')
title('CMAES: blue-avg red-best','interpreter','latex')
ylim([0 1])

% legend({'Average fitnesses','Best fitnesses'},'interpreter','latex','location','southeast')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])

cd CMAES
cd(runfilename)
saveas(gcf,[runfilename '_fitevoltrial.fig'],'fig')
saveas(gcf,[runfilename '_fitevoltrial.eps'],'epsc')
saveas(gcf,[runfilename '_fitevoltrial.png'],'png')
cd ..
cd ..


%% best shape
parameters{4}=1;
alpha=griddecoder_cont(res.CMAES.bestgenome);
parameters{4}=0;
title(['CMAES best shape: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

cd CMAES
cd(runfilename)
saveas(gcf,[runfilename '_bestshape.fig'],'fig')
saveas(gcf,[runfilename '_bestshape.eps'],'epsc')
saveas(gcf,[runfilename '_bestshape.png'],'png')
cd ..
cd ..

%% flatshape
parameters{4}=1;
alpha=griddecoder_cont(ones(size(res.CMAES.bestgenome)));
parameters{4}=0;
title(['Flat layer: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

cd CMAES
cd(runfilename)
saveas(gcf,[runfilename '_fullshape.fig'],'fig')
saveas(gcf,[runfilename '_fullshape.eps'],'epsc')
saveas(gcf,[runfilename '_fullshape.png'],'png')
cd ..
cd ..



%% Cutoff study
% % meshplotter(parameters{1},bestgenome(:,end))
% % temp=popgenome(:,1);
% for cutoff=0:10:100
%     temp=xmin;
%     [sorted,index]=sort(temp);
%     % temp=popgenome(:,10);
%     % cutoff=0;
%     temp(index(1:cutoff))=0;
%     temp(index(cutoff+1:100))=1;
%     % meshplotter(parameters{1},popgenome(:,10))
%     % meshplotter(parameters{1},temp)
%
%     parameters{4}=1;
%     objfunc(temp)
%     cd CMAES
%     saveas(gcf,['CMAES_cutoff_' num2str(cutoff) '.fig'],'fig')
%     saveas(gcf,['CMAES_cutoff_' num2str(cutoff) '.eps'],'epsc')
%     saveas(gcf,['CMAES_cutoff_' num2str(cutoff) '.png'],'png')
%     cd ..
% end

%% best shape from each trial without discretising 

figure(15)
for trial=1:n_trials
    clf

discgenome = res.CMAES.bestgenometrial(:,trial);

alpha=griddecoder_cont(discgenome);
meshplotter(mesh,discgenome);

title(['CMAES best shape before rounding: trial ' num2str(trial) ' fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)
% pause


cd CMAES
cd(runfilename)
if trial==1
mkdir bestshapetrial
end
cd bestshapetrial
saveas(gcf,[runfilename '_bestshapetrial' num2str(trial) '.fig'],'fig')
saveas(gcf,[runfilename '_bestshapetrial' num2str(trial) '.eps'],'epsc')
saveas(gcf,[runfilename '_bestshapetrial' num2str(trial) '.png'],'png')
cd ..
cd ..
cd ..
end

%% best shape from each trial after discretising 
% 
% figure(15)
% for trial=1:n_trials
%     clf
% % parameters{4}=1;
% % alpha=griddecoder_cont(res.CMAES.bestgenome);
% discgenome = zeros(size(res.CMAES.bestgenometrial(:,trial)));
% discgenome(res.CMAES.bestgenometrial(:,trial)>0.5)=1;
% % discgenome = floor(2*res.CMAES.bestgenometrial(:,trial));
% alpha=griddecoder_comb(discgenome);
% meshplotter(mesh,discgenome);
% % parameters{4}=0;
% title(['CMAES best shape after rounding: trial ' num2str(trial) ' fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)
% % pause
% 
% 
% cd CMAES
% cd(runfilename)
% if trial==1
% mkdir bestshapetrial
% end
% cd bestshapetrial
% saveas(gcf,[runfilename '_bestshapetrial' num2str(trial) '_disc.fig'],'fig')
% saveas(gcf,[runfilename '_bestshapetrial' num2str(trial) '_disc.eps'],'epsc')
% saveas(gcf,[runfilename '_bestshapetrial' num2str(trial) '_disc.png'],'png')
% cd ..
% cd ..
% cd ..
% end

