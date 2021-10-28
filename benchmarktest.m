% test benchmark
global BM videotoggle
videotoggle=0;

for BM=[7]
clear
close all
global BM;
disp(BM)
prob=benchmarks(BM);
disp(prob.problem)
AFSO_input_from_benchmark;

% MAIN
runfname=['benchmark_' num2str(BM) '_long']; % GA DE CMAES
budget=20; n_trials=1;
N = mesh.N;           % number of objective variables/problem dimension
npop = 4+floor(3*log(N));  % population size, offspring number

%     % TEST
%     runfname=['benchmark_' num2str(BM) '_test']; % GA DE CMAES
%     npop=1; budget=100; n_trials=7;
%
%     % for SIMP & SIMP_wof
%     runfname=['benchmark_' num2str(BM) '_long'];
%     npop=32; budget=4096; n_trials=31;
%
%     % for SIMPr
%     runfname=['benchmark_' num2str(BM) '_long'];
%     npop=1; budget=200; n_trials=31;
%
%     runfname=['benchmark_' num2str(BM) '_'];
%     npop=32; budget=2048; n_trials=31;
%
%     % for TEST RUNS
%     runfname=['benchmark_' num2str(BM) '_xi_0.5'];
%     npop=5; budget=100; n_trials=1;
    
    %% Algorithms
    run_HC
%     run_GA
%     run_tabu
%     run_CH
%     run_CH_D
%     run_CH_R
%     run_DE
%     run_DE_D
%     run_CMAES
%     run_CMAES_D
%     run_CMAES_P
%     budget=200;
%     run_SIMP
%     run_SIMPr
%     run_SIMP_wof
%     run_HC
%     run_PF
%     run_PUFIELD
%     run_PUFIELD2
%     run_PUFIELD3
%     run_MOHC
%     run_MOHC_TCH
%     run_MOCMAESD
%     run_MOCMAESDf2
%     run_MOCMAES
%     run_MOTABU
%     run_MOGA
%     run_MOSIMP

%     run_SIMP_HC
    
    
    
end
return

%% plotting average fitness evolution
% Algorithms={'GA','DE','CMAES','CH'};
figure
Algorithms={'GA','CMAES_D','CH','CH_D'};
markers={'b-','r--','k-.','m','g','k'};
hold on
for alg=1:length(Algorithms)
    cd(Algorithms{alg})
    runfilename=[runfname Algorithms{alg}];
    cd(runfilename)
    load(runfilename,'res')
    algores=getfield(res,Algorithms{alg});
    BESTFITEVOLTRIAL{alg}=algores.bestfitevol;
    FEVALS{alg}=algores.fevalsgen;
    plot(FEVALS{alg},BESTFITEVOLTRIAL{alg},markers{alg},'linewidth',2,'displayname',[Algorithms{alg} ' best fitness evol']);
    cd ..
    cd ..
end
set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])
xlabel('Function evaluations','interpreter','latex')
ylabel('Fitness','interpreter','latex')
title('Avg Best Fitness Evolution across trials','interpreter','latex')
ylim auto
legend('location','southeast','interpreter','none')
grid on
box on

cd benchmarks
mkdir(runfname)
cd(runfname)
saveas(gcf,[runfname '_bestfitevol.fig'],'fig')
saveas(gcf,[runfname '_bestfitevol.eps'],'epsc')
saveas(gcf,[runfname '_bestfitevol.png'],'png')
cd ..
cd ..


%% plotting whisker plot of the distribution of best fitnesses from each algorithm
figure
X=[];G=[];
for alg=1:length(Algorithms)
    cd(Algorithms{alg})
    runfilename=[runfname Algorithms{alg}];
    cd(runfilename)
    load(runfilename,'res')
    algores=getfield(res,Algorithms{alg});
    BESTFITTRIAL{alg}(1,:)=algores.bestfitnesstrial;
    FEVALS{alg}=algores.fevalsgen;
    X=[X BESTFITTRIAL{alg}];
    G=[G alg*ones(size(BESTFITTRIAL{alg}))];
    
    labels{alg}=Algorithms{alg};
    
    cd ..
    cd ..
end

boxplot(X,G,'labels',labels);


title('GA','interpreter','latex')
% ylim([0 1])
ylim auto
set(findobj(gca,'type','line'),'linew',1.5)
% set(gca,'TickLabelInterpreter','tex','LineWidth',1.5,'Fontsize',16)
set(gca,'LineWidth',1.5,'Fontsize',16)
set(gcf ,'position',[10   100   800 600],'Color',[1 1 1])

xlabel('Algorithms','interpreter','latex')
ylabel('Fitness','interpreter','latex')
title('Distribution of Best fitness trials','interpreter','latex')

% legend('location','southeast','interpreter','latex')
grid on
box on

cd benchmarks
mkdir(runfname)
cd(runfname)
saveas(gcf,[runfname '_bestfitdist.fig'],'fig')
saveas(gcf,[runfname '_bestfitdist.eps'],'epsc')
saveas(gcf,[runfname '_bestfitdist.png'],'png')
cd ..
cd ..







