%% construction heuristic for the shape optimisation problem

runfilename=[runfname 'CH'];
% load(runfilename);

tic
begintime=toc;

%% Memoization
objfunc=@griddecoder_comb;
% objfunc=memoize(@griddecoder_comb);
% objfunc.CacheSize=100;
% stat=stats(objfunc);

%% Initiating savefile
cd CH
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




%% construction heuristic controls

V0=mesh.N; %  mesh.N:no volume constraint
N=mesh.N;
totalfeval=mesh.N*(mesh.N+1)/2 - (mesh.N-V0)*(mesh.N-V0+1)^2/2;
disp(['Total fevals: ' num2str(totalfeval)]) 
remporous=V0; % no of porous materials remaining to be filled.
filled=zeros(mesh.N,1); % initially nothing is filled
remaining=find(filled==0);fevalcounter=0;

numbertofill=ceil(V0*(2*N-V0)/(2*budget-(2*N-V0)));


% numbertofill=1;

disp(['Filling ' num2str(numbertofill) ' in each construction step']);
estfevals=sum(length(mesh.domain):-numbertofill:0);
disp(['Estimated fevals : ' num2str(estfevals)])
g=0;
genmax=ceil(V0/numbertofill);
cutoff=zeros(1,genmax);

while remporous>0
    g=g+1;
    maxfitness=0;
    ifitness=zeros(length(remaining),1);
    for i=1:length(remaining) % fill in i location

        genome=zeros(mesh.N,1); % reset genome
        genome(filled==1)=1; % fill porous in the filled places
        %         remaining=find(filled==0);
        genome(remaining(i))=1;
        fitness=objfunc(genome);
        ifitness(i)=fitness;
        fevalcounter=fevalcounter+1;
        
        %         %% GIF WRITER
        %         filename=['Constructionanimation' num2str(frequencies) '.gif' ];
        %         title([ 'Construction heuristic sweep: ' num2str(frequencies)  'Hz (Fitness=' num2str(fitness) ')'],'FontSize',15)
        %
        %         % place this outside for loop
        %         % place the following inside for loop after plotting each image
        %         % for i=1:...
        %         frame = getframe(gcf);
        %         im = frame2im(frame);
        %         [imind,cm] = rgb2ind(im,256);
        %         % Write to the GIF File
        %         if fevalcounter == 1  % the beginning of the loop
        %             imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
        %         else
        %             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
        %         end
        %         %% GIF WRITER end
        if fitness>maxfitness
            maxi=i;
            maxfitness=fitness;
        end
        if mod(i,10)==0
            %% calculating percent complete and time remaining
            curriteration=fevalcounter;
            totiteration=estfevals;
            currenttime=toc;percent_complete=curriteration/(totiteration)*100;
            remaining_time_estimator(currenttime,curriteration,totiteration,percent_complete)
        end
    end
    res.CH.fevalsgen(1,g)=fevalcounter;

    
    
    [ascendfitness,ascendindex]=sort(ifitness);
    if length(ascendindex)<numbertofill
        tofill=ascendindex;
    else
        tofill=ascendindex(end:-1:end-numbertofill+1);
    end
    [maxfitnessitr,maxloc]=max(ifitness);
    

    res.CH.bestgenevol(:,g)=filled;
    res.CH.bestgenevol(maxloc,g)=1;
        res.CH.bestfitevol(1,g)=maxfitnessitr;
    if g>1 && res.CH.bestfitevol(1,g-1)>res.CH.bestfitevol(1,g)
        res.CH.bestfitevol(1,g)=res.CH.bestfitevol(1,g-1);
        res.CH.bestgenevol(:,g)=res.CH.bestgenevol(:,g-1);
    end
        
        
        
    filled(remaining(tofill))=1; % filling in the maximum absorption location
    
    cutoff(1,g)=sum(filled);
    res.CH.bestgenomes(:,g)=filled;
    res.CH.bestfitnesses(1,g)=objfunc(filled);
    
    
    remaining=find(filled==0);
    remporous=V0-sum(filled);
    disp(['Filled ' num2str(sum(filled)) ' out of ' num2str(length(filled))])
    
    %% saving data file
    cd CH
    cd(runfilename)
    save(runfilename)
    disp('Results are saved in ')
    disp(runfilename)
    cd ..
    cd ..
    %% calculating percent complete and time remaining
    curriteration=fevalcounter;
    totiteration=totalfeval;
    currenttime=toc;percent_complete=curriteration/(totiteration)*100;
    remaining_time_estimator(currenttime,curriteration,totiteration,percent_complete)
end

[bestfitness,maxindex]=max(res.CH.bestfitnesses);
bestgenome=res.CH.bestgenomes(:,maxindex);

res.CH.bestgenome=bestgenome;
res.CH.bestfitness=bestfitness;

currenttime=toc;
disp(['Time for completion = ' num2str(currenttime) ' s'])

%% saving data

cd CH
cd(runfilename)
save(runfilename)
disp('Results are saved in ')
disp(runfilename)
cd ..
cd ..


%% Plotting average fitness evolution

figure(17)
box on

% plot(res.CH.fevalsgen,res.CH.avgfitevol,'b-','linewidth',2,'displayname','Average fitness in population')

hold on
plot(res.CH.fevalsgen,res.CH.bestfitevol,'r--','linewidth',2,'displayname','Best fitness in population')

xlabel('Function evaluations','interpreter','latex')
ylabel('Fitness ','interpreter','latex')
title('CH','interpreter','latex')
ylim([0 1])
legend('location','southeast','interpreter','latex')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf ,'position',[10   100   800 600],'Color',[1 1 1])


cd CH
cd(runfilename)
saveas(gcf,[runfilename '_fitevol.fig'],'fig')
saveas(gcf,[runfilename '_fitevol.eps'],'epsc')
saveas(gcf,[runfilename '_fitevol.png'],'png')
cd ..
cd ..


%% Best shape

figure
meshplotter(mesh,res.CH.bestgenome)
title(['CH best shape : Fitness = ' num2str(res.CH.bestfitness,3)],'interpreter','latex','Fontsize',16)
set(gcf ,'position',[10   100   800 600],'Color',[1 1 1])
cd CH
cd(runfilename)
saveas(gcf,[runfilename '_bestshape.fig'],'fig')
saveas(gcf,[runfilename '_bestshape.eps'],'epsc')
saveas(gcf,[runfilename '_bestshape.png'],'png')
cd ..
cd ..




%% Loop through each generation fitnesses for best shapes (animation)

for i=1:length(res.CH.bestfitnesses)
    clf
    meshplotter(mesh,res.CH.bestgenomes(:,i))
    weight=sum(res.CH.bestgenomes(:,i));
    title(['CH best shape : Fitness = ' num2str(res.CH.bestfitnesses(i),3) ' no of 1s: ' num2str(weight) ],'interpreter','latex','Fontsize',16)


    keyboard
end




