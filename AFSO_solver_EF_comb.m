function [alpha,SAC_2mic]=AFSO_solver_EF_comb(mesh,frequences,varargin)
% [alpha,SACAbhi]=AFSO_solver_EF_comb(mesh,frequencies,mat,plotToggle,order,displayToggle,PressureFieldToggle)
% % This is a solver which computes the alpha given the mesh and frequencies
% assumes Equivalent fluid model (no poroelasticity)
% alpha= average absorption coefficient over all frequencies
% SAC= Sound absorption coefficient at each frequency in the form of array

% mesh is a structure containing all mesh info
% frequences is an array of frequencies for which alpha is requested
% mat is a structure containing material property array and material model
% plotToggle=1 shows absorption curve after solving using 2mic & imp method
% FEorder indicates whether the mesh info is 1st or 2nd order or more
% displayToggle=1 prints additional information in console while code runs
% PressureFieldToggle=1 animates pressure field after solving

% This code has contributions from
% Abhilash , Savvas (FE solver)
% Luc Jaouen, Fabien Chevillote (models semi, airprophr)
% Vivek (fast sparse assembly, i/o, bc, 2mic method, derivatives)
% Compiled and extended by: Vivek (vivek23991@gmail.com)
% absorption curves verified with alphacell 
% Date: March 27, 2020. 
% disp('BH initiated')
%% extracting data
if nargin>=3
    mat=varargin{1};
else
    error('input incomplete, material property missing')
end
if nargin>=4
    plotToggle=varargin{2};
else
    plotToggle=0;
end

if nargin>=5
    order=varargin{3};
else
    order=1;
end
if nargin>=6
    displayToggle=varargin{4};
else
    displayToggle=0;
end
if nargin>=7
    PressureFieldToggle=varargin{5};
else
    PressureFieldToggle=0;
end

%% Hard coded constants
 % SIMP PENALTY 
xiMin = 1E-03 ; % default interpolation value for air since Xi can't be 0.
% xiMin should be close to zero.
tolelemfrac=0.1; % tolerance fraction of element size for boundary node extraction
thickness = 1; %  thickness of imp tube outsid plane of view
SparseIndex=1; % whether sparse matrix should be used
FA=1; % whether faster assembly procedure should be used.
elembw2mic=1; % the number of air element between 2 microphones
%% --------------------calculated values----------------------
omega = 2*pi*frequences ;

%% Getting impedance tube dimensions automatically

h_sample = 2*57e-3; % [m] sample thickness %  this value is unused in this code
xLim = [min(mesh.coords(:,1)), max(mesh.coords(:,1))] ; % [xmin xmax] of the impedance tube
yLim = [min(mesh.coords(:,2)), max(mesh.coords(:,2))] ; % [ymin ymax] of the impedance tube
chamberHeight = abs(yLim(2) - yLim(1));

%% Meshing geometry
tol=min(tolelemfrac*mesh.D/mesh.NELEMxD,tolelemfrac*mesh.d/mesh.NELEMy); 
meshType = mesh.type;%'TRI' ;

% Element type Xi = 0 for air and Xi=1 for poroelastic
Xi = ones(length(mesh.connect(:,1)),1) ;

% Creating air elements
Xi(mesh.matType==1) = xiMin ;
%% Creating boundary conditions
% NOTE: For equivalent fluid model, the side walls being symmetric and
% having rigid wall are the same boundary condition. It is implicit and
% doesnt have to be specified in the FE code
% Pressure boundaries
presBCDofs = [];
pres = 1;
% Front edge (left side) excited by acoustic wave of unit amplitude
leftNodes = find(abs(mesh.coords(:,1)-xLim(1))<tol) ;
presBCNodes = leftNodes ;
presBCDofs = [presBCDofs ; presBCNodes ] ;
% presBCDofs = presBCDofs + 2*size(mesh.coords,1); % accounts for ordering of Dofs

% Collecting all Dirichlet boundary DOFs and free DOFs
constrDofs = sort(presBCDofs) ;
freeDofs = setdiff(1:length(mesh.coords(:,1)),constrDofs) ;
%% ------------ Macroscopic parameters-------------%
switch order
    case 1
        if(strcmp(meshType,'QUAD'))
            ElementsType = 10; % 2D poroelastic quadrilateral element
            nGP = 4; % number of Gauss Points for numerical integration
            nodeFE = (order+1)^2; % number of nodes per quadrilateral finite-element
        elseif(strcmp(meshType,'TRI'))
            ElementsType = 8; % 2D poroelastic triangle element
            nGP = 3 ;
            nodeFE = 3 ;
        end
    case 2
        if(strcmp(meshType,'QUAD'))
            ElementsType = 20; % 2D poroelastic quadrilateral element
            nGP = (order+1)^2; % number of Gauss Points for numerical integration
            nodeFE = (order+1)^2; % number of nodes per quadrilateral finite-element
        elseif(strcmp(meshType,'TRI'))
            ElementsType = 82; % 2D poroelastic quadrilateral element
            nGP = 7; % number of Gauss Points for numerical integration
            nodeFE = 6; % number of nodes per quadrilateral finite-element
        end
end

E=mat.prop(10);%1.6E05 ;
nu = mat.prop(11);%0.44;
%% Run from here if you dont wanna change the mesh

% frequences = 200:200:4000;
omega = 2*pi*frequences ;

% mat = [0.99 , 1.96E-04, 9.8E-05, 1E04, 1.01, 1.399,1.25, 4.75E-09, 8.0, 1.6E05, 0.44, 0.1];
% [eta_s,rho_tilde,~,gamma_tilde,~,~,rho_eq,K_eq,~,zair ] =matPropsBiotTemp2 (mat.prop, frequences,mat.model,h_sample/2) ;
[eta_s,rho_tilde,~,gamma_tilde,~,~,rho_eq,K_eq,~,zair ] =matPropsBiot_vivek (mat.prop, frequences,mat.model,h_sample/2) ;

% Creating analogous material parameters for air
temperature_celsius=20;   %[25]
pression_Pa=101325;          %[101325]
humidite_relative=45;     %[80]
[rho_Air,c,~,gamma,~,~]=propairhr(temperature_celsius,pression_Pa,humidite_relative);
K_Air=gamma*pression_Pa;
kAir = omega / c; % wave number of air

rhoInterp = rho_Air+xiMin*(rho_eq-rho_Air) ;
KInterp = K_Air + xiMin*(K_eq-K_Air) ;

%% ------------FE preliminaries-------------%
dofFE = nodeFE ; % number of disp dofs per element
nTotNodes = length(mesh.coords(:,1)); % total number of nodes
nTotdofs = nTotNodes; % total number of disp nodes
nTotElem = length(mesh.connect(:,1)); % total number of elements
thickness = 1; % unit thickness

Hp= sparse(nTotNodes,nTotNodes); Ha= sparse(nTotNodes,nTotNodes);
Qp = sparse(nTotNodes,nTotNodes); Qa = sparse(nTotNodes,nTotNodes);
    

%% ------------Assembly-------------%

% Loop through each finite-element
for i=1:nTotElem
    if mod(i,100)==0 && displayToggle==1
        disp('Elements assembled:')
        disp(i)
    end
    %Coords of element
    XYZ(1,:)=mesh.coords(mesh.connect{i},1);
    XYZ(2,:)=mesh.coords(mesh.connect{i},2);
    
    % Helem and Qelem computed without material parameters
    Helem = Get_Element_H_BiotAcoustics(ElementsType,XYZ,thickness,nGP,nodeFE);
    Qelem = Get_Element_Q_BiotAcoustics(ElementsType,XYZ,thickness,nGP,nodeFE);
    
    % Assembling into global matrices
    
    if(Xi(i)==xiMin)
        %         [ Ha, Qa ] = StiffnessMatrixFEM_HelmholtzAcoustics( i,Ha,Helem,Qa,Qelem,mesh) ;
        % modified code by vivek begins
        sctrA=cell2mat(mesh.connect(i,:))';
        Ha(sctrA,sctrA)=Ha(sctrA,sctrA)+Helem;
        Qa(sctrA,sctrA)=Qa(sctrA,sctrA)+Qelem;
        % modified code by vivek ends
    elseif(Xi(i)==1)
        %         [ Hp, Qp ] = StiffnessMatrixFEM_HelmholtzAcoustics( i,Hp,Helem,Qp,Qelem,mesh) ;
        % modified code by vivek begins
        sctrA=cell2mat(mesh.connect(i,:))';
        Hp(sctrA,sctrA)=Hp(sctrA,sctrA)+Helem;
        Qp(sctrA,sctrA)=Qp(sctrA,sctrA)+Qelem;
        % modified code by vivek ends
    end
    %     if mod(i,10)==0
    %         disp(['Assembling elements: ' num2str(100*i/nTotElem,4) '% Complete..'])
    %     end
    
end
%%
%------------Incporating frequency dependent parts-------------%
% i = sqrt(-1) ;
SACAbhi=zeros(length(frequences),1);
Zsn=zeros(length(frequences),1);
for iFreq = 1:length(frequences)
    %     disp(['Frequency under computation: ' num2str(frequences(iFreq),5) 'Hz'])
    F = sparse(zeros(nTotNodes,1));
    %     HTilde = Hp/rho_eq(iFreq) + Ha/rhoInterp(iFreq);
    %     QTilde = Qp/K_eq(iFreq) + Qa/KInterp(iFreq);
    HTilde = (1/rho_eq(iFreq))*sparse(Hp) + (1/rhoInterp(iFreq))*sparse(Ha);
    QTilde = 1/K_eq(iFreq)*sparse(Qp) + 1/KInterp(iFreq)*sparse(Qa);
    
    StiffTilde = HTilde - omega(iFreq)^2*QTilde;
    
    % Applying Dirichlet-type BCs
    sol = zeros(nTotNodes,1) ; % 3 Dofs per node
    sol(presBCDofs) = pres ; % enforcing pressures on left boundary
    F = F -StiffTilde*sol ;
    sol(freeDofs) = full((StiffTilde(freeDofs,freeDofs)) \ (F(freeDofs))) ; % solving only for unconstrained cases
    %% plotting pressure field
%         PressureFieldToggle=0;
    if iFreq==1 && PressureFieldToggle==1
        keyboard
        %% surface plot
%         figure
%         for theta=0:(2*pi/60):4*pi
%             clf
%             
%             hold on
%             for iElem=1:nTotElem
%                 nodes=mesh.connect{iElem};
%                 x=mesh.coords(nodes,1);
%                 y=mesh.coords(nodes,2);
%                 z=real(sol(nodes)*exp(-1i*theta));
%                 
%                 %                 fill3(x,y,z)
%                 %                 patch('Faces',nodes,'Vertices',mesh.coords(nodes),'FaceVertexCData',z,'FaceColor','interp');
%                 patch(x,y,z,'FaceAlpha',.5)
%                 if mesh.matType(iElem)==2
%                 %    patch('Faces',[1 2 3 4],'Vertices',mesh.coords(nodes,:),'FaceColor','red','FaceAlpha',.5)
%                     patch('Faces',[1 2 3 4],'Vertices',mesh.coords(nodes,:),'EdgeColor','black','FaceColor','none','LineWidth',2,'FaceColor','black','FaceAlpha',.5)
%                 end
%                 
%                 
% 
%             end
%             pause(0.001)
%         end
        %% line plot
                    nodeshori=(mesh.NELEMxD+mesh.NELEMxL+1);
            
            midrow=round((nTotNodes/nodeshori)/2);
            
            nodes=((midrow-1)*nodeshori+1):((midrow)*nodeshori);
            ylimit=[-max(abs(sol)) max(abs(sol))];
            
        for theta=0:(2*pi/60):4*pi
            clf

            plot(mesh.coords(nodes,1),real(sol(nodes)*exp(1i*theta)))
                        ylim(ylimit)
                        title(['Pressure field along the centreline' ''])
            pause(0.1)

        end
    end     % END OF PRESSURE FIELD PLOTS

    
    
    
    
    %% two microphone method to compute absorption
    x1=mesh.L/3; % x coordinate of first microphone
    x2=mesh.L/3+elembw2mic*mesh.L/mesh.NELEMxL; % x coordinate of second microphone (the very next element)
    y=mesh.d/2; % y coordinate of both the microphones, here its the centre of impedance tube
%     dx=x1-x2 is the distance between the two microphones
    [~,x1_node]=min((mesh.coords(:,1)-x1).^2+(mesh.coords(:,2)-y).^2);
    [~,x2_node]=min((mesh.coords(:,1)-x2).^2+(mesh.coords(:,2)-y).^2);
    
    if x1_node==x2_node
        error('the two microphones are placed in the same node, please check')
    end
    
    Px1_dof=x1_node; % dof corresponding to pressure at microphone 1
    Px2_dof=x2_node; % dof corresponding to pressure at microphone 2
    P1=sol(Px1_dof); % complex fourier pressure amplitude at microphone 1
    P2=sol(Px2_dof); % complex fourier pressure amplitude at microphone 2
    k=2*pi*frequences(iFreq)/c; % wave number in air
%     Rn=(P1*exp(-1i*k*x2)-P2*exp(-1i*k*x1) ) / (P1*exp(1i*k*x2)-P2*exp(1i*k*x1)) ; % -B/A LKKK
    Rn=(P1*exp(-1i*k*x2)-P2*exp(-1i*k*x1) ) / (-P1*exp(1i*k*x2)+P2*exp(1i*k*x1)); % B/A Luc
    SAC_2mic(iFreq) = 1-(abs(Rn))^2;
    
    %
    %     % to get normal fluid displacements
    Fg = StiffTilde * sol/omega(iFreq).^2 ;
    Un = sum(Fg(presBCDofs) )/ chamberHeight;
    Zsn = -pres / (Un*1i*omega(iFreq) * zair)  ;
    
    SAC_imp(iFreq)=1-abs((Zsn -1)/(Zsn +1))^2;
    
    
end
%% Final alpha
% alpha=mean(SAC_imp);
alpha=mean(SAC_2mic);

%% plotting
if(plotToggle==1)
    % Visualization symmetries of displacement
    
    figure(23)
    set(gca,'FontSize',16)
    hold on
    grid on
    
    plot(frequences,SAC_imp,'b--x','Linewidth',2,'Displayname','Impedance method') ; % FEM soln
    plot(frequences,SAC_2mic,'r--o','Linewidth',2,'Displayname','Two microphone method') ; % FEM soln
    xlabel('Frequency (Hz)')
    ylabel('\alpha - absorption coefficient')
    
    title('Equivalent fluid model')
    set(gca,'Ylim',[0 1])
    legend
end

%% sanity check: prints warning if absorption values are negative
if sum(SAC_imp<0)>0 || sum(SAC_2mic<0)>0
    warning('Some absorption values are negative!')
    disp('One cause could be less number of air layers in front')
    %     keyboard
end
if mean(abs(SAC_2mic-SAC_imp))*100>5
    warning(['Average percentage error bw 2mic and p/u  ' num2str(mean(abs(SAC_2mic-SAC_imp))*100) ' % '])
% keep an eye on this difference as this can be very large at low
% frequencies
end
