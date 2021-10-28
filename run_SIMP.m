%%  SIMP from benchmark

runfilename=[runfname 'SIMP'];
% load(runfilename);
mesh.penalty=3;

%% Memoization
objfunc=@griddecoder_cont; %(maximisation)
% objfunc=memoize(@griddecoder_cont);
% objfunc.CacheSize=100;
% stat=stats(objfunc);

%% Initiating savefile
cd SIMP
mkdir(runfilename)
cd(runfilename)
if isfile([runfilename '.mat'])
    warning('A save file with the same name already exists. Rename or create backup')
%     keyboard
else
    save(runfilename)
end
cd ..
cd ..

%% Restart if previous run exists
checkPreviousRuns

%% SIMP controls
nelx=mesh.NELEMxD;
nely=mesh.NELEMy;
volfrac=1;%mesh.VF;
penal=mesh.penalty;
rmin=2;
ft=2;
genmax=budget;
npop=1;
move = 0.2;

%% Random number seed
rng(123)

res.SIMP.trial(n_trials).bestgenomes=zeros(genomelength,genmax);
res.SIMP.trial(n_trials).bestfitnesses=zeros(1,genmax);
% res.SIMP.trial(n_trials).popevol=zeros(genomelength,npop,genmax);
res.SIMP.trial(n_trials).fitevol=zeros(genmax,1);
res.SIMP.avgfitevoltrial=zeros(genmax,n_trials);
res.SIMP.bestfitevoltrial=zeros(genmax,n_trials);



tic
begintime=toc;

for trial=1:n_trials
    disp(['Trial : ' num2str(trial) ' started'])

    %% Trial based Random number seed
    disp('Using random number seed based on trial number')
    rngseed=123+(trial-1)*1000000;
    rng(rngseed)
    %% ACTUAL ALGORITHM
    bestfitnesses=zeros(1,genmax);
    bestgenomes=zeros(genomelength,genmax);
    % PREPARE FILTER
    iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (i1-1)*nely+j1;
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                    e2 = (i2-1)*nely+j2;
                    k = k+1;
                    iH(k) = e1;
                    jH(k) = e2;
                    sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                end
            end
        end
    end
    H = sparse(iH,jH,sH); % turns out to be identity matrix for some reason
    Hs = sum(H,2); % the row sums of H matrix as a vector
    %% INITIALIZE ITERATION
    if trial==1
        %         xinit = repmat(volfrac,nely,nelx); % 0.5 in all
        xinit = repmat(0.5,nely,nelx); % 0.5 in all
    else
        
        xinit = rand(nely,nelx); % random solution
%         xinit=xinit./mean(xinit)*volfrac; % normalizing based on volume fraction
%         xinit(xinit>1)=1;
%         xinit(xinit<0)=0;
    end
    
    
    x=xinit;
    xPhys = x;
%     g = 0;
    change = 1;
    %% START ITERATION
%     x_itr=zeros(nely,nelx,45);
%     while change > 0.01 && g<genmax
%         g = g + 1;
    for g=1:genmax
        
        %% FE-ANALYSIS % evaluation of the objective function presteps
        %     sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        %     K = sparse(iK,jK,sK); K = (K+K')/2;
        %     U(freedofs) = K(freedofs,freedofs)\F(freedofs);
        %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        %     ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
        %     c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
        %     dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce; % derivative of obj func
        % those elements which are in the design domain
        currgenome=zeros(mesh.NELEMy*mesh.NELEMxD,1);
        for irow=1:mesh.NELEMy
            row=xPhys(end-irow+1,:);
            currgenome(mesh.NELEMxD*(irow-1)+[1:mesh.NELEMxD])=row';
        end
        mesh.matType(mesh.domain)=currgenome+1; % setting 0 in genome to 1 (air) and 1 in genome to 2(porous)
        
        %         x_itr(:,:,g)=xPhys;
        
        %% plotting each gen
        %         figure(1)
        %         clf
        %         meshplotter(mesh,currgenome)
        %         title(['MMA iteration: ' num2str(g) ' SAC: calculating...'  ' Vol: ' num2str(mean(xPhys(:))) ] )
        
        
        
        %% solving
        
        [alpha,~,da_dXi]=AFSO_solver_BH_cont_adjoint(mesh,parameters{2},parameters{7},parameters{3},parameters{5},parameters{6},parameters{9});
        
        c=(1-alpha);%+sum(xPhys.*(1-xPhys));
        
        bestgenomes(:,g)=currgenome;
        bestfitnesses(:,g)=alpha;
        
        deriv_mat=zeros(size(x));
        for irow=1:mesh.NELEMy
            row=da_dXi(mesh.NELEMxD*(irow-1)+[1:mesh.NELEMxD]);
            deriv_mat(mesh.NELEMy-irow+1,:)=row;
        end
        dc=-deriv_mat;%+(1-2*xPhys);
        
        
        dv = ones(nely,nelx);
        %% FILTERING/MODIFICATION OF SENSITIVITIES
        if ft == 0
            % no filtering
            
        elseif ft==1
            % first filtering technique
            dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
            
        elseif ft == 2
            % second filtering technique
            dc(:) = H*(dc(:)./Hs);
            dv(:) = H*(dv(:)./Hs);
        end
        %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        l1 = 0; l2 = 1e9; % move = 0.2;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            %         xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));% this is yielding complex numbers
            xnew = max(0,max(x-move,min(1,min(x+move,real(x.*sqrt(-dc./dv/lmid))))));% this is yielding complex numbers
            if ft==0
                xPhys=xnew;
            elseif ft == 1
                xPhys = xnew;
            elseif ft == 2
                xPhys(:) = (H*xnew(:))./Hs;
            end
            %% plotting temporary mesh
            %         figure(2)
            %          colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
            
            %% % % %
            if sum(xPhys(:)) > volfrac*nelx*nely
                l1 = lmid;
            else
                l2 = lmid;
            end
        end
        change = max(abs(xnew(:)-x(:)));
        x = xnew;
        %% PRINT RESULTS
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',g,c, ...
            mean(xPhys(:)),change);
%         colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
        
        
        
        %% PLOT DENSITIES
        %     colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
        % figure(1)
        % clf
        % meshplotter(mesh,genome)
        %         title(['MMA itr: ' num2str(g) ' SAC: ' num2str(alpha,'%11.4f') ' Vol: ' num2str(mean(xPhys(:))) ] )
        % keyboard
        
        %     cd mma
        %     saveas(1,[runfilename 'itr' num2str(loop) 'shape.fig'],'fig')
        %
        %     saveas(1,[runfilename 'itr' num2str(loop) 'shape.eps'],'epsc')
        %     saveas(1,[runfilename 'itr' num2str(loop) 'shape.png'],'png')
        %     % save([runfilename 'itr' num2str(loop) 'genome.mat'],'genome','alpha','da_dXi')
        %     cd ..
        % place the following inside for loop after plotting each image
        % place the following inside for loop after plotting each image
        % for i=1:...
        %% GIF WRITER
        cd SIMP
        cd(runfilename)
        if trial==1 || trial==2 || trial==3 || trial==4
            writegif=0;
        else
            writegif=0;
        end
        if writegif==1
            meshplotter(mesh,currgenome)
            title(['SIMP trial: ' num2str(trial) ' itr: ' num2str(g) ' fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

            % place this outside for loop
            giffilename=[runfilename 'trial' num2str(trial) 'gif.gif'];
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if g == 1  % change this variable name and value to the beginning of the loop
                imwrite(imind,cm,giffilename,'gif', 'Loopcount',inf,'DelayTime',.5);
            else
                imwrite(imind,cm,giffilename,'gif','WriteMode','append','DelayTime',.5);
            end
        end
        cd ..
        cd ..
        [~,maxfitnessgeneration]=max(round(bestfitnesses,4));
        if change<0.01 || maxfitnessgeneration < g-10 % last best was before 5 generations
            disp('TERMINATION DUE TO CONVERGENCE')
            for gg=g+1:size(bestfitnesses,2)
                bestgenomes(:,gg)=bestgenomes(:,maxfitnessgeneration);
            end
            bestfitnesses(:,g+1:end)=bestfitnesses(maxfitnessgeneration);
            break
        end
    end
    %% saving data in structure
    
    fitevol=bestfitnesses;
%      res.SIMP.trial(trial).bestgenomes=bestgenomes; % bestgenome in generation
%         res.SIMP.trial(trial).bestfitnesses=bestfitness; % best fitnesses in generation
    [res.SIMP.bestfitnesstrial(1,trial),maxindex]=max(bestfitnesses);
    res.SIMP.bestgenometrial(:,trial)=bestgenomes(:,maxindex);
    
    %     res.SIMP.trial(trial).popevol=popevol;
    res.SIMP.trial(trial).fitevol=fitevol;
    
    res.SIMP.avgfitevoltrial(:,trial)=res.SIMP.trial(trial).fitevol;%mean(fitevol,2); % each column is avg fitness across generation in each trial
    res.SIMP.avgfitevol=mean(res.SIMP.avgfitevoltrial(:,1:trial),2); % avg across trials of avg fitness evolution
    
    res.SIMP.bestfitevoltrial(:,trial)=res.SIMP.trial(trial).fitevol;%max(res.SIMP.trial(trial).fitevol,[],2);
    res.SIMP.bestfitevol=mean(res.SIMP.bestfitevoltrial(:,1:trial),2);
    
    [res.SIMP.trial(trial).bestfitness,maxindex]=max(bestfitnesses);
    res.SIMP.trial(trial).bestgenome=bestgenomes(:,maxindex);
    

    res.SIMP.bestgenomedtrial(:,trial)=round(res.SIMP.bestgenometrial(:,trial));
    res.SIMP.bestfitnessdtrial(1,trial)=griddecoder_comb(res.SIMP.bestgenomedtrial(:,trial));

    %% saving mat file after finishing the trial
    cd SIMP
    cd(runfilename)
    save(runfilename)
    cd ..
    cd ..
    disp(['Trial ' num2str(trial) ' finished'])
    


    
    %% progress display
    OPS_runProgressTracker
end
%------------------------------

%% post proc

[res.SIMP.bestfitness,besttrial]=max(res.SIMP.bestfitnesstrial);
res.SIMP.bestgenome=res.SIMP.bestgenometrial(:,besttrial);
res.SIMP.fevalsgen=1:1:genmax*1;
disp('Hurray! SIMP completed.')

%% discretizing best genomes and shapes

res.SIMP.bestgenomedtrial=round(res.SIMP.bestgenometrial);
for trial=1:n_trials
res.SIMP.bestfitnessdtrial(1,trial)=griddecoder_comb(res.SIMP.bestgenomedtrial(:,trial));
disp(['trial ' num2str(trial) ' discretized! '])
end

[res.SIMP.bestfitnessd,besttriald]=max(res.SIMP.bestfitnessdtrial);
res.SIMP.bestgenomed=res.SIMP.bestgenomedtrial(:,besttriald);

%% Save

cd SIMP
cd(runfilename)
save(runfilename)
cd ..
cd ..

endtime=toc;
disp(['Total time elapsed: ' num2str(endtime)])

%% Plotting average fitness evolution

figure(17)
box on

plot(res.SIMP.fevalsgen,res.SIMP.avgfitevol,'b-','linewidth',2,'displayname','Average fitness in population')

hold on
plot(res.SIMP.fevalsgen,res.SIMP.bestfitevol,'r--','linewidth',2,'displayname','Best fitness in population')

xlabel('Function evaluations','interpreter','latex')
ylabel('Average fitness across trials','interpreter','latex')
title('SIMP','interpreter','latex')
ylim([0 1])
legend('location','southeast','interpreter','latex')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])


cd SIMP
cd(runfilename)
saveas(gcf,[runfilename '_fitevol.fig'],'fig')
saveas(gcf,[runfilename '_fitevol.eps'],'epsc')
saveas(gcf,[runfilename '_fitevol.png'],'png')
cd ..
cd ..

%% plotting fitness evolution for all trials
figure(18)

% plot(res.SIMP.avgfitevoltrial,'b','linewidth',1,'displayname','Average fitness in population')

% hold on

plot(res.SIMP.bestfitevoltrial,'linewidth',1,'displayname','Best fitness in population')

xlabel('Generation number','interpreter','latex')
ylabel('Trialwise fitnesses','interpreter','latex')
title('SIMP: all trials','interpreter','latex')
ylim([0 1])

% legend({'Average fitnesses','Best fitnesses'},'interpreter','latex','location','southeast')

set(gca,'TickLabelInterpreter','latex','LineWidth',1.5,'Fontsize',16)
set(gcf,'Color',[1 1 1])

cd SIMP
cd(runfilename)
saveas(gcf,[runfilename '_fitevoltrial.fig'],'fig')
saveas(gcf,[runfilename '_fitevoltrial.eps'],'epsc')
saveas(gcf,[runfilename '_fitevoltrial.png'],'png')
cd ..
cd ..


%% best shape
parameters{4}=1;
trial=1;
meshplotterlite(mesh,res.SIMP.trial(trial).bestgenome)
alpha=res.SIMP.trial(trial).bestfitness;

title(['SIMP best shape: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

title(['SIMP trial ' num2str(trial) ' best shape: fitness = ' num2str(alpha,3)],'interpreter',  'latex','fontsize',16)

cd SIMP
cd(runfilename)
% saveas(gcf,[runfilename '_bestshape.fig'],'fig')
% saveas(gcf,[runfilename '_bestshape.eps'],'epsc')
% saveas(gcf,[runfilename '_bestshape.png'],'png')

saveas(gcf,[runfilename ' trial ' num2str(trial) '_bestshape.fig'],'fig')
saveas(gcf,[runfilename ' trial ' num2str(trial) '_bestshape.eps'],'epsc')
saveas(gcf,[runfilename ' trial ' num2str(trial) '_bestshape.png'],'png')
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
%     cd SIMP
%     saveas(gcf,['SIMP_cutoff_' num2str(cutoff) '.fig'],'fig')
%     saveas(gcf,['SIMP_cutoff_' num2str(cutoff) '.eps'],'epsc')
%     saveas(gcf,['SIMP_cutoff_' num2str(cutoff) '.png'],'png')
%     cd ..
% end






