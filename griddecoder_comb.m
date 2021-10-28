function [fitness,SACAbhi]=griddecoder_comb(genome,varargin)
% passing function for the mesh decoder in the shape optimisation project
% This version assumes only one half of the domain and uses boundary
% conditions for a symmetric FE problem by considering no slip in the impedance
% tube centreline
% The syntax for this funcion is
% [fitness,abscurve]=griddecoder1(genome,parameters)
% genome is a vector containing according number of elements for
% determining the choice of porous or air - '0' means air and '1' means
% porous
% griddecoder1 can solve both symmetric and full models
if sum((genome>0.001).*(genome<0.999))
   error('genome sent to grid decoder_comb should have only either 0s or 1s')
end
global parameters 

%% extraction of variables from parameters
mesh=parameters{1};
frequencies=parameters{2};
plotToggle=parameters{3};
meshPlotToggle=parameters{4};
FEorder=parameters{5};
displayToggle=parameters{6};
mat=parameters{7};
skeletonmodel=parameters{8};
pressureFieldToggle=parameters{9};

%% editing the mesh to represent the element types according to the genome

switch mesh.bc
    case 'full'
%         mesh.domain=mesh.matType==2; % those elements which are in the design domain
        symmgenome=[genome; genome]; % generating a symmetric genome for the upper half
        % note that symmgenome is not yet symmetric. The latter half is made
        % symmetric in the following for loop.
        
        for i=1:mesh.NELEMy/2 % for each row in the mesh
            row=genome(mesh.NELEMxD*(i-1)+1:mesh.NELEMxD*i); % copying ith row in genome
            % finding which row in the symmgenome this should be copied to
            j=mesh.NELEMy-i+1; % the symmetric row number is tot no of rows -i+1
            symmgenome(mesh.NELEMxD*(j-1)+1:mesh.NELEMxD*j)=row;
        end
        mesh.matType(mesh.domain)=symmgenome+1; % setting 0 in genome to 1 (air) and 1 in genome to 2(porous)
        
    case 'symmetric'
%         mesh.domain=find(mesh.matType==2); % those elements which are in the design domain
        mesh.matType(mesh.domain)=genome+1; % setting 0 in genome to 1 (air) and 1 in genome to 2(porous)
        
end



%% solving by passing to AFSO_solver
switch skeletonmodel
    case 'BiotHelmholtz'
  [fitness,SACAbhi]=AFSO_solver_BH_comb(mesh,frequencies,mat,plotToggle,FEorder,displayToggle,pressureFieldToggle);

    case 'EquivalentFluid'
    [fitness,SACAbhi]=AFSO_solver_EF_comb(mesh,frequencies,mat,plotToggle,FEorder,displayToggle,pressureFieldToggle); % this function evaluates the problem and gives the fitness
    otherwise
    error('Incorrect model name')
end



%% printing output
printToggle=0;
if printToggle
    global runfilename
    global currentcount
    global AFSODATA
    printfilename=[runfilename '.txt'];
    problemID=[ (mesh.NELEMxD) (mesh.NELEMy) (mesh.D) (mesh.d) frequencies(1) frequencies(end) length(frequencies)  strcmp(mesh.bc,'full') strcmp(mesh.sidebc,'fixed')  ];
    % num2str(problemID')
    % probleminstance=[ num2str(mesh.NELEMxD) ' x ' num2str(mesh.NELEMy) ' grid ' num2str(mesh.D) ' m x' num2str(mesh.d) ' m ' num2str(frequencies(1)) ' to ' num2str(frequencies(end)) ' Hz '  num2str(length(frequencies)) ' steps ' mesh.bc ' problem ' mesh.sidebc ' sidewall' ];
    
    probleminstance=[ num2str(problemID(1)) ' x '  num2str(problemID(2)) ' grid ' num2str(problemID(3)) ' m x' num2str(problemID(4)) ' m ' num2str(problemID(5)) ' to ' num2str(problemID(6)) ' Hz '  num2str(problemID(7)) ' steps '  num2str(problemID(8)) ' fullproblem '  num2str(problemID(9)) ' fixedsidewall' ];
    materialcard=[ mat.model ' ' num2str(mat.prop(1)) ' '  num2str(mat.prop(2)) ' '  num2str(mat.prop(3)) ' '  num2str(mat.prop(4)) ' '  num2str(mat.prop(5)) ' '  num2str(mat.prop(6)) ' '  num2str(mat.prop(7)) ' '  num2str(mat.prop(8)) ' '  num2str(mat.prop(9)) ' '  num2str(mat.prop(10)) ' ' num2str(mat.prop(11)) ' ' num2str(mat.prop(12))  ];
    
    solution=[];
    for ii=1:length(genome)
        solution=[solution num2str(genome(ii))];
    end
    solution=[solution ' ' num2str(fitness)];
    % fprintf(['\n' probleminstance '\n'])
    % fprintf([materialcard '\n'])
    
    %=length(AFSODATA.fitness);
    AFSODATA.genome{currentcount+1}=genome;
    AFSODATA.fitness{currentcount+1}=fitness;
    AFSODATA.problemID{currentcount+1}=problemID;
    AFSODATA.materialcard{currentcount+1}=materialcard;
    currentcount=currentcount+1;
    cd landscape
    save([runfilename '.mat'],'AFSODATA','-append')
    fid1=fopen(printfilename,'a+');
    fprintf(fid1,'%s\r\n',solution); % print the line with end line
    fclose(fid1);
    cd ..
end

%% plot the mesh
if meshPlotToggle==1
    figure(10)
    clf
    hold on  
    meshplotter(parameters{1},genome)
    title(['Mesh fitness = ' num2str(fitness)],'Interpreter','Latex','FontSize',23)
    pause(0.1)
end







