%% deconstruction heuristic for the shape optimisation problem

runfilename=[runfname 'CH_D'];
% load(runfilename);

tic
begintime=toc;

%% Memoization
% objfunc=@griddecoder_comb;
objfunc=memoize(@griddecoder_comb);
objfunc.CacheSize=100;
stat=stats(objfunc); 

%% Initiating savefile
cd CH_D
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




%% deconstruction heuristic controls

V0=mesh.N; %  mesh.N:no volume constraint
N=mesh.N;
totalfeval=mesh.N*(mesh.N+1)/2 - (mesh.N-V0)*(mesh.N-V0+1)^2/2;
disp(['Total fevals with 1 notofill: ' num2str(totalfeval)]) 
remair=V0; % no of porous materials remaining to be filled.
filled=ones(mesh.N,1); % initially nothing is filled
remaining=find(filled==1);fevalcounter=0;

notofill=ceil(V0*(2*N-V0)/(2*budget-(2*N-V0)));
% notofill=1;
totalfevalcorr= sum(mesh.N:-notofill:0);

disp(['Total fevals with ' num2str(notofill) ' notofill: ' num2str(totalfevalcorr)]) 
g=0;

genmax=ceil(V0/notofill);
cutoff=zeros(1,genmax);
% res.CH_D.bestgenomes=zeros(genomelength,genmax);
% res.CH_D.bestgenevol=zeros(genomelength,genmax);
% res.CH_D.fevalsgen=zeros(1,genmax);

while remair>0
    g=g+1;
    maxfitness=0;
    ifitness=zeros(length(remaining),1);
    for i=1:length(remaining) % fill in i location

        genome=ones(mesh.N,1); % reset genome
        genome(filled==0)=0; % fill air in the air places
        %         remaining=find(filled==0);
        genome(remaining(i))=0;
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
            totiteration=totalfevalcorr;
            percent_complete=curriteration/(totiteration)*100;
            
            currenttime=toc;
            ETR=currenttime/curriteration*(totiteration)-currenttime;
            ETRday=floor(ETR/(24*3600));
            ETRhr=floor(mod(ETR,24*3600)/3600);
            ETRmin=floor(mod(ETR,3600)/60);
            ETRsec=mod(ETR,60);
            disp(['---------- Remaining time estimate------------'])
            disp(['Avg time per feval: ' num2str(currenttime/curriteration)])
            disp(['Est. time rem: ' num2str(ETRday) ' D ' num2str(ETRhr) ' H ' num2str(ETRmin) ' M ' num2str(ETRsec) ' S'])
            disp(['Percent complete : ' num2str(percent_complete) ' % '])
            disp(['----------------------------------------------'])
        end
    end
    res.CH_D.fevalsgen(1,g)=fevalcounter;

    
    
    [ascendfitness,ascendindex]=sort(ifitness);
    tofill=ascendindex(end:-1:max(end-notofill+1,1));
    [maxfitnessitr,maxloc]=max(ifitness);
    

%     res.CH_D.bestgenevol(:,g)=filled;
%     res.CH_D.bestgenevol(maxloc,g)=0;
        res.CH_D.bestfitevol(1,g)=maxfitnessitr;
    if g>1 && res.CH_D.bestfitevol(1,g-1)>res.CH_D.bestfitevol(1,g)
        res.CH_D.bestfitevol(1,g)=res.CH_D.bestfitevol(1,g-1);
%         res.CH_D.bestgenevol(:,g)=res.CH_D.bestgenevol(:,g-1);
    end
        
        
        
    filled(remaining(tofill))=0; % removing porous in the maximum absorption location
    
    cutoff(1,g)=sum(filled);
    res.CH_D.bestgenomes(:,g)=filled;
    res.CH_D.bestfitnesses(1,g)=objfunc(filled);
    
    
    remaining=find(filled==1);
    remair=sum(filled);
    disp(['Porous remaining ' num2str(sum(filled)) ' out of ' num2str(length(filled))])
    
    %% saving data file
    cd CH_D
    cd(runfilename)
    save(runfilename)

    disp(runfilename)
    cd ..
    cd ..
    
    
    
    %% calculating percent complete and time remaining
    curriteration=fevalcounter;
    totiteration=totalfevalcorr;
    percent_complete=curriteration/(totiteration)*100;
    
    currenttime=toc;
    ETR=currenttime/curriteration*(totiteration)-currenttime;
    ETRday=floor(ETR/(24*3600));
    ETRhr=floor(mod(ETR,24*3600)/3600);
    ETRmin=floor(mod(ETR,3600)/60);
    ETRsec=mod(ETR,60);
    disp(['---------- Remaining time estimate------------'])
    disp(['Avg time per feval: ' num2str(currenttime/curriteration)])
    disp(['Est. time rem: ' num2str(ETRday) ' D ' num2str(ETRhr) ' H ' num2str(ETRmin) ' M ' num2str(ETRsec) ' S'])
    disp(['Percent complete : ' num2str(percent_complete) ' % '])
    disp(['----------------------------------------------'])
end

[bestfitness,maxindex]=max(res.CH_D.bestfitnesses);
bestgenome=res.CH_D.bestgenomes(:,maxindex);

res.CH_D.bestgenome=bestgenome;
res.CH_D.bestfitness=bestfitness;

currenttime=toc;
disp(['Time for completion = ' num2str(currenttime) ' s'])

%% saving data

cd CH_D
cd(runfilename)
save(runfilename)
disp('Results are saved in ')
disp(runfilename)
cd ..
cd ..


%% Plotting average fitness evolution

figure(17)
box on

% plot(res.CH_D.fevalsgen,res.CH_D.avgfitevol,'b-','linewidth',2,'displayname','Average fitness in population')

hold on
plot(res.CH_D.fevalsgen,res.CH_D.bestfitevol,'r--','linewidth',2,'displayname','Best fitness in population')

xlabel('Function evaluations','interpreter','latex')
ylabel('Fitness ','interpreter','latex')
title('CH\_D','interpreter','latex')
ylim([0 1])
legend('location','southeast','interpreter','latex')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])


cd CH_D
cd(runfilename)
saveas(gcf,[runfilename '_fitevol.fig'],'fig')
saveas(gcf,[runfilename '_fitevol.eps'],'epsc')
saveas(gcf,[runfilename '_fitevol.png'],'png')
cd ..
cd ..


%% Best shape

figure
meshplotter(mesh,res.CH_D.bestgenome)
title(['CH\_D best shape : Fitness = ' num2str(res.CH_D.bestfitness,3)],'interpreter','latex','Fontsize',16)
cd CH_D
cd(runfilename)
saveas(gcf,[runfilename '_bestshape.fig'],'fig')
saveas(gcf,[runfilename '_bestshape.eps'],'epsc')
saveas(gcf,[runfilename '_bestshape.png'],'png')
cd ..
cd ..




%% Loop through each generation fitnesses for best shapes
% 
% for i=1:length(res.CH_D.bestfitnesses)
%     clf
%     meshplotter(mesh,res.CH_D.bestgenomes(:,i))
%     title(['CH\_D best shape : Fitness = ' num2str(res.CH_D.bestfitnesses(i),3)],'interpreter','latex','Fontsize',16)
% 
% 
%     pause(0.1)
% end




