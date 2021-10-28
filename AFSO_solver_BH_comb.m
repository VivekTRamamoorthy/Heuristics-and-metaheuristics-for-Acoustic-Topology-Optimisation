function [alpha,SAC_2mic]=AFSO_solver_BH_comb(mesh,frequences,varargin)
% [alpha,SACAbhi]=AFSO_solver_BH_comb(mesh,frequencies,mat,plotToggle,order,displayToggle,PressureFieldToggle)
% % This is a solver which computes the alpha given the mesh and frequencies
% assumes Biot Helmholtz model (poroelastic model)
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
% p=1; % SIMP PENALTY 
xiMin = 1E-03 ; % default interpolation value for air since Xi can't be 0.
% xiMin should be close to zero.
toleranceelemfrac=0.1; % tolerance fraction of element size for boundary node extraction
thickness = 1; %  thickness of imp tube outsid plane of view
SparseIndex=1; % whether sparse matrix should be used
FA=1; % whether faster assembly procedure should be used.
elembw2mic=1; % the number of air element between 2 microphones
%% calculated values
omega = 2*pi*frequences ;

%% Getting impedance tube dimensions automatically

h_sample = 2*57e-3; % [m] sample thickness %  this value is unused in this code
xLim = [min(mesh.coords(:,1)), max(mesh.coords(:,1))] ; % [xmin xmax] of the impedance tube
yLim = [min(mesh.coords(:,2)), max(mesh.coords(:,2))] ; % [ymin ymax] of the impedance tube
chamberHeight = abs(yLim(2) - yLim(1));

%% Meshing geometry
meshType = mesh.type;% 'TRI' or 'QUAD' ;
% tolerance for extracting nodes for applying the boundary conditions
tol=min(toleranceelemfrac*mesh.D/mesh.NELEMxD,0.1*mesh.d/mesh.NELEMy);
% tol =toleranceelemfrac * element width or height whichever is minimum
% Element type Xi = 0 for air and Xi=1 for poroelastic
Xi = ones(length(mesh.connect(:,1)),1) ;

% Creating air elements
Xi(mesh.matType==1) = xiMin ;
%% CREATING BOUNDARY CONDITIONS
dispBCDofs = [];

% Right edge  fully bonded
rightNodes = find(abs(mesh.coords(:,1)-xLim(2))<tol) ;
dispBC_rearDofs = [dispBCDofs ; sort([2*rightNodes-1;2*rightNodes]) ] ;

% BOTTOM EDGE BOUNDARY CONDITION
bottomNodes = find(abs(mesh.coords(:,2)-yLim(1))<tol) ;   % Bottom
% Choice of the boundary condition for the side walls of impedance tube
if isfield(mesh,'sidebc')
    sidebc=mesh.sidebc;
else
    sidebc='symmetric'; disp('Side walls assumed symmetric')
    % sidebc='fixed'; disp('Side walls assumed fixed')
end
switch sidebc
    case 'symmetric'     % only fixing the uy degrees of freedom
        dispBC_bottomDofs = [dispBCDofs ; 2*bottomNodes ] ;
    case 'fixed' % fixing both the ux and uy degrees of freedom
        dispBC_bottomDofs = [dispBCDofs ; sort([2*bottomNodes-1;2*bottomNodes])] ;
    otherwise
        error('check the sidebc. It should be symmetric or fixed')
end

% TOP EDGE BOUNDARY CONDITION
topNodes = find(abs(mesh.coords(:,2)-yLim(2)) <tol );   % Top
% Top edge BC varies basec on whether the FE model is for only bottom half and top is symmetric
if strcmp(mesh.bc,'symmetric') % if FE model is only for bottom half
    % only fixing the ux degree of freedom
    dispBC_topDofs = [dispBCDofs ; 2*topNodes ] ;
elseif strcmp(mesh.bc,'full') % if FE model is for full impedance tube
    switch sidebc
        case 'symmetric'
            % only fixing the ux degrees of freedom
            dispBC_topDofs = [dispBCDofs ; 2*topNodes ] ;
        case 'fixed'
            % fixing both the ux and uy degrees of freedom
            dispBC_topDofs = [dispBCDofs ; sort([2*topNodes-1;2*topNodes])] ;
        otherwise
            error('check the sidebc. It should be symmetric or fixed')
    end
end
    
% collecting all left nodes: Left nodes are load nodes and no bc is applied
leftNodes = find(abs(mesh.coords(:,1)-xLim(1))<tol) ;

% Collecting all constrained displacement DOFs
dispBCDofs = sort([ dispBC_rearDofs ; dispBC_bottomDofs ; dispBC_topDofs ] );

% Pressure 0 boundaries
pres0BCDofs = [];
pres0 = 0;
% No nodes have pressure =0 in this problem so we leave it empty
pres0BCNodes = [] ;
pres0BCDofs = [pres0BCDofs ; pres0BCNodes ] ;
pres0BCDofs = pres0BCDofs + 2*size(mesh.coords,1); % accounts for ordering of Dofs



% Pressure boundaries
presBCDofs = [];
pres = 1;
% Front edge (left side) excited by acoustic wave of unit amplitude

presBCNodes = leftNodes ;
presBCDofs = [presBCDofs ; presBCNodes ] ;
presBCDofs = presBCDofs + 2*size(mesh.coords,1); % accounts for ordering of Dofs

% Collecting all Dirichlet boundary DOFs and free DOFs
constrDofs = sort([dispBCDofs ; presBCDofs; pres0BCDofs]) ;
freeDofs = setdiff(1:3*length(mesh.coords(:,1)),constrDofs) ;

%%
%------------ Macroscopic parameters-------------%
switch order
    case 1
        if(strcmp(meshType,'QUAD'))
%             ElementsType = 10; % 2D poroelastic quadrilateral element
            nGP = 4; % number of Gauss Points for numerical integration
            nodeFE = (order+1)^2; % number of nodes per quadrilateral finite-element
            loc=[-0.577350269189626 , -0.577350269189626;
                0.577350269189626 , -0.577350269189626;
                0.577350269189626 ,  0.577350269189626;
                -0.577350269189626 ,  0.577350269189626];
            wt=[1;1;1;1];
        elseif(strcmp(meshType,'TRI'))
%             ElementsType = 8; % 2D poroelastic triangle element
            nGP = 3 ;
            nodeFE = 3 ;
        end
    case 2
        if(strcmp(meshType,'QUAD'))
            ElementsType = 20; % 2D poroelastic quadrilateral element
            nGP = (order+1)^2; % number of Gauss Points for numerical integration
            nodeFE = (order+1)^2; % number of nodes per quadrilateral finite-element
            loc=[-0.774596669241483 , -0.774596669241483;
                0.774596669241483 , -0.774596669241483;
                0.774596669241483 ,  0.774596669241483
                -0.774596669241483 ,  0.774596669241483;
                0.000000000000000 , -0.774596669241483;
                0.774596669241483 ,  0.000000000000000;
                0.000000000000000 ,  0.774596669241483;
                -0.774596669241483 ,  0.000000000000000;
                0.000000000000000 ,  0.000000000000000];
            
            r1wt(1,1) = 0.555555555555556;
            r1wt(2,1) = 0.555555555555556;
            r1wt(3,1) = 0.888888888888889;
            wt=[r1wt(1,1)*r1wt(1,1);
                r1wt(1,1)*r1wt(1,1);
                r1wt(1,1)*r1wt(1,1);
                r1wt(1,1)*r1wt(1,1);
                r1wt(1,1)*r1wt(3,1);
                r1wt(1,1)*r1wt(3,1);
                r1wt(1,1)*r1wt(3,1);
                r1wt(1,1)*r1wt(3,1);
                r1wt(3,1)*r1wt(3,1)];
        elseif(strcmp(meshType,'TRI'))
            ElementsType = 82; % 2D poroelastic quadrilateral element
            nGP = 7; % number of Gauss Points for numerical integration
            nodeFE = 6; % number of nodes per quadrilateral finite-element
        end
end

E=mat.prop(10);%1.6E05 ;
nu = mat.prop(11);%0.44;

%% Run from here if you dont wanna change the mesh


[eta_s,rho_tilde,~,gamma_tilde,~,~,rho_eq,K_eq,~,zair,c ] =matPropsBiot_vivek (mat.prop, frequences,mat.model,h_sample/2) ;
E_por = E*(1+1i*eta_s) ;
%%  Creating analogous material parameters for air
temperature_celsius=20;   %[25]
pression_Pa=101325;          %[101325]
humidite_relative=45;     %[80]
[rho_Air,~,~,gamma,~,~]=propairhr(temperature_celsius,pression_Pa,humidite_relative);
K_Air=gamma*pression_Pa;
E_Air = 0 ;
nu_Air = 0 ;

% Interpolating the material parameters
rhoInterp = rho_Air+xiMin*(rho_eq-rho_Air) ;
KInterp = K_Air + xiMin*(K_eq-K_Air) ;
EInterp = E_Air + Xi*(E_por-E_Air) ;
nuInterp = nu_Air + Xi*(nu-nu_Air) ;
rho_s_tildeInterp = xiMin.*rho_tilde ;
gamma_tildeInterp = xiMin.*gamma_tilde ;
DeInterp{length(Xi)}=[];

for ii=1:length(Xi)
    coeff = EInterp(ii)/((1+nuInterp(ii))*(1-2*nuInterp(ii))) ;
    DeInterp{ii} = coeff*[1-nuInterp(ii) , nuInterp(ii) , 0 ; nuInterp(ii) , 1-nuInterp(ii) ,  0 ; 0 ,0 , (1-2*nuInterp(ii))/2] ;
end
%%
%------------FE preliminaries-------------%
dofFE = 2*nodeFE ; % number of disp dofs per element
nTotNodes = length(mesh.coords(:,1)); % total number of nodes
nTotdofs = 2 * nTotNodes; % total number of disp nodes
nTotElem = length(mesh.connect(:,1)); % total number of elements
if(SparseIndex==0)
    Kp = zeros(nTotdofs); Ka = zeros(nTotdofs);
    Mp = zeros(nTotdofs); Ma = zeros(nTotdofs);
    Hp= zeros(nTotNodes); Ha= zeros(nTotNodes);
    Qp = zeros(nTotNodes); Qa = zeros(nTotNodes);
    Cp = zeros(nTotdofs, nTotNodes); Ca = zeros(nTotdofs, nTotNodes);
    Fs = zeros(nTotdofs,1);
    Fp = zeros(nTotNodes,1);
elseif(SparseIndex==1)
    Kp = sparse(nTotdofs,nTotdofs); Ka = sparse(nTotdofs,nTotdofs);
    Mp = sparse(nTotdofs,nTotdofs); Ma = sparse(nTotdofs,nTotdofs);
    Hp= sparse(nTotNodes,nTotNodes); Ha= sparse(nTotNodes,nTotNodes);
    Qp = sparse(nTotNodes,nTotNodes); Qa = sparse(nTotNodes,nTotNodes);
    Cp = sparse(nTotdofs, nTotNodes); Ca = sparse(nTotdofs, nTotNodes);
    Fs = sparse(nTotdofs,1);
    Fp = sparse(nTotNodes,1);
end
%% ------------Assembly-------------%
% disp('Assembling...')
% tic
if FA==1
    %% faster assembly
    % initializing global vectors for rows and columns of the sparse global
    % matrices
    G_ROW_A=zeros(nTotElem*(length(mesh.connect{1}))^2,1);
    G_ROW_B=zeros(nTotElem*(2*length(mesh.connect{1}))^2,1);
    G_ROW_BA=zeros(nTotElem*2*(length(mesh.connect{1}))^2,1);
    
    
    G_COL_A=zeros(nTotElem*(length(mesh.connect{1}))^2,1);
    G_COL_B=zeros(nTotElem*(2*length(mesh.connect{1}))^2,1);
    G_COL_BA=zeros(nTotElem*2*(length(mesh.connect{1}))^2,1);
    
    % initializing global vectors for the elements of global matrices
    
    Ka1=zeros(size(G_ROW_B));
    Ma1=zeros(size(G_ROW_B));
    Qa1=zeros(size(G_ROW_A));
    Ha1=zeros(size(G_ROW_A));
    Ca1=zeros(size(G_ROW_BA));
    
    Kp1=zeros(size(G_ROW_B));
    Mp1=zeros(size(G_ROW_B));
    Qp1=zeros(size(G_ROW_A));
    Hp1=zeros(size(G_ROW_A));
    Cp1=zeros(size(G_ROW_BA));
    
end
% Loop through each finite-element
for i=1:nTotElem
    if mod(i,100)==0 && displayToggle==1
        disp('Elements assembled:')
        disp(i)
    end
    %Coords of element
    XYZ(1,:)=mesh.coords(mesh.connect{i},1);
    XYZ(2,:)=mesh.coords(mesh.connect{i},2);
    
    Kelem = zeros(dofFE , dofFE) ; Melem = zeros(dofFE , dofFE) ; Helem = zeros(nodeFE , nodeFE) ; Qelem = zeros(nodeFE , nodeFE) ; Celem = zeros(dofFE, nodeFE);
    
    for j=1:nGP
        gp=loc(j,:);
        xi=gp(1); eta=gp(2); zeta = 0 ;
        %Shape functions and Derivatives of shape functions
        switch order
            case 1
                ShapeFunction=1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta) ;
                    (1+xi)*(1+eta);(1-xi)*(1+eta)];
                DerivativesOfShapeFunction=1/4*[-(1-eta), 1-eta, 1+eta, -(1+eta) ;
                    -(1-xi) ,-(1+xi),1+xi,    1-xi];
            case 2
                N1 = xi*(1-xi)*eta*(1-eta)/4; N3 = -xi*(1+xi)*eta*(1-eta)/4; N5 = xi*(1+xi)*eta*(1+eta)/4; N7 = -xi*(1-xi)*eta*(1+eta)/4;
                N2 = -(1+xi)*(1-xi)*eta*(1-eta)/2; N4 = xi*(1+xi)*(1-eta)*(1+eta)/2; N6 = (1+xi)*(1-xi)*eta*(1+eta)/2; N8 = -xi*(1-xi)*(1-eta)*(1+eta)/2; N9 = (1-xi^2)*(1-eta^2);
                ShapeFunction=[N1;N2;N3;N4;N5;N6;N7;N8;N9];
                
                dN1deta=xi*(1-xi)*(1-2*eta)/4; dN3deta=-xi*(1+xi)*(1-2*eta)/4; dN5deta=xi*(1+xi)*(1+2*eta)/4; dN7deta=-xi*(1-xi)*(1+2*eta)/4; dN2deta=-(1+xi)*(1-xi)*(1-2*eta)/2;
                dN4deta=xi*(1+xi)*(-2*eta)/2; dN6deta=(1+xi)*(1-xi)*(1+2*eta)/2; dN8deta=-xi*(1-xi)*(-2*eta)/2; dN9deta=(1+xi)*(1-xi)*(-2*eta);
                dN1dxi=(1-2*xi)*eta*(1-eta)/4; dN3dxi=-(1+2*xi)*eta*(1-eta)/4; dN5dxi=(1+2*xi)*eta*(1+eta)/4; dN7dxi=-(1-2*xi)*eta*(1+eta)/4; dN2dxi=2*xi*eta*(1-eta)/2;
                dN4dxi=(1+2*xi)*(1-eta)*(1+eta)/2; dN6dxi=-2*xi*eta*(1+eta)/2; dN8dxi=-(1-2*xi)*(1-eta)*(1+eta)/2; dN9dxi=-2*xi*(1+eta)*(1-eta);
                
                dNdxi=[dN1dxi;dN2dxi;dN3dxi;dN4dxi;dN5dxi;dN6dxi;dN7dxi;dN8dxi;dN9dxi];
                dNdeta=[dN1deta;dN2deta;dN3deta;dN4deta;dN5deta;dN6deta;dN7deta;dN8deta;dN9deta];
                DerivativesOfShapeFunction=[dNdxi';dNdeta'];
        end
        %Calculate the Jacobian
        Jacob = DerivativesOfShapeFunction*XYZ' ;
        GradNpmicro = inv(Jacob)*DerivativesOfShapeFunction;
        Npmicro    = ShapeFunction';
        N = zeros(2,dofFE);
        N(1,1:2:end)   = ShapeFunction';
        N(2,2:2:end)   = ShapeFunction';
        %Calculate the Volume of element
        dVolume=det(Jacob)*wt(j,1)*thickness;
        %Calculate B1 Matrix
        B1 = zeros(3,4) ; detJacob=Jacob(2,2)*Jacob(1,1)-Jacob(1,2)*Jacob(2,1);
        B1= [Jacob(2,2) , Jacob(1,2) , 0 , 0 ; 0 , 0 , -Jacob(2,1) , Jacob(1,1) ; -Jacob(2,1) , Jacob(1,1) , Jacob(2,2) , -Jacob(1,2)]/detJacob ;
        %Calculate B2 matrix
        B2=zeros(4,2*size(DerivativesOfShapeFunction,2));
        for k=1:size(DerivativesOfShapeFunction,2)
            B2(1,2*k-1)=DerivativesOfShapeFunction(1,k);
            B2(2,2*k-1)=DerivativesOfShapeFunction(2,k);
            B2(3,2*k)=DerivativesOfShapeFunction(1,k);
            B2(4,2*k)=DerivativesOfShapeFunction(2,k);
            k=k+1;
        end
        %Calculate B matrix
        [ B ] = B1*B2;
        %Stiffness matrix of RVE element
        Kelem=Kelem+B'*DeInterp{i}*B*dVolume;
        Melem = Melem + N'*N*dVolume;
        Helem = Helem + GradNpmicro'*GradNpmicro*dVolume;
        Qelem = Qelem + Npmicro'*Npmicro*dVolume;
        Celem = Celem + N'*GradNpmicro*dVolume;
    end
    
    %Pressure Field
    sctrA = mesh.connect{i} ;
    %Displacement Field
    sctrB=zeros(1,2*length(sctrA));
    sctrB(1:2:end) = 2*sctrA-1;
    sctrB(2:2:end) = 2*sctrA;
    %% faster assembly
    if FA==1
        ROW_A=sctrA'*ones(size(sctrA)); % forming the row number matrix
        ROW_B=sctrB'*ones(size(sctrB));
        ROW_BA=sctrB'*ones(size(sctrA));
        
        COL_A=ones(size(sctrA'))*sctrA; % forming the column number matrix
        COL_B=ones(size(sctrB'))*sctrB;
        COL_BA=ones(size(sctrB'))*sctrA;
        
        totelemA=length(ROW_A(:));
        totelemB=length(ROW_B(:));
        totelemBA=length(ROW_BA(:));
        
        RANGE_A=(totelemA*(i-1)+1) : (totelemA*i); % location to copy in vector
        RANGE_B=(totelemB*(i-1)+1) : (totelemB*i);
        RANGE_BA=(totelemBA*(i-1)+1) : (totelemBA*i);
        
        G_ROW_A(RANGE_A)=ROW_A(:); % global row vector
        G_ROW_B(RANGE_B)=ROW_B(:);
        G_ROW_BA(RANGE_BA)=ROW_BA(:);
        
        G_COL_A(RANGE_A)=COL_A(:); % global column vector
        G_COL_B(RANGE_B)=COL_B(:);
        G_COL_BA(RANGE_BA)=COL_BA(:);
    end
    %%
    
    if(Xi(i)==xiMin)
        if FA==0
            Ka(sctrB,sctrB) = Ka(sctrB,sctrB) + Kelem;
            Ma(sctrB,sctrB) = Ma(sctrB,sctrB) + Melem;
            Ha(sctrA,sctrA) = Ha(sctrA,sctrA) + Helem;
            Qa(sctrA,sctrA) = Qa(sctrA,sctrA) + Qelem;
            Ca(sctrB,sctrA) = Ca(sctrB,sctrA) + Celem;
        else
            
            %% fast assembly
            Ka1(RANGE_B)=Kelem(:);
            Ma1(RANGE_B) =  Melem(:);
            Ha1(RANGE_A) =  Helem(:);
            Qa1(RANGE_A) = Qelem(:);
            Ca1(RANGE_BA) = Celem(:);
            
        end
        
        
        
    elseif(Xi(i)==1)
        if FA==0
            Kp(sctrB,sctrB) = Kp(sctrB,sctrB) + Kelem;
            Mp(sctrB,sctrB) = Mp(sctrB,sctrB) + Melem;
            Hp(sctrA,sctrA) = Hp(sctrA,sctrA) + Helem;
            Qp(sctrA,sctrA) = Qp(sctrA,sctrA) + Qelem;
            Cp(sctrB,sctrA) = Cp(sctrB,sctrA) + Celem;
        else
            %% fast assembly
            Kp1(RANGE_B)=Kelem(:);
            Mp1(RANGE_B) =  Melem(:);
            Hp1(RANGE_A) =  Helem(:);
            Qp1(RANGE_A) = Qelem(:);
            Cp1(RANGE_BA) = Celem(:);
            
        end
    end
end
% end of assembly
if FA==1
    %%
    Ka=sparse(G_ROW_B,G_COL_B,Ka1);
    Ma=sparse(G_ROW_B,G_COL_B,Ma1);
    Ha=sparse(G_ROW_A,G_COL_A,Ha1);
    Qa=sparse(G_ROW_A,G_COL_A,Qa1);
    Ca=sparse(G_ROW_BA,G_COL_BA,Ca1);
    
    Kp=sparse(G_ROW_B,G_COL_B,Kp1);
    Mp=sparse(G_ROW_B,G_COL_B,Mp1);
    Hp=sparse(G_ROW_A,G_COL_A,Hp1);
    Qp=sparse(G_ROW_A,G_COL_A,Qp1);
    Cp=sparse(G_ROW_BA,G_COL_BA,Cp1);
    
    
end
% toc




%%------------Incporating frequency dependent parts-------------%
%%
KTilde = Kp + Ka  ;

% SACAbhi=zeros(length(frequences),1);%Zsn=zeros(length(frequences),1);

for iFreq = 1:length(frequences)
    if displayToggle==1
        disp([' Computing for frequency : ' num2str(frequences(iFreq))]);
    end
    F = [Fs;Fp];
    MTilde = Mp * rho_tilde(iFreq) + Ma * rho_s_tildeInterp(iFreq) ;
    HTilde = Hp/rho_eq(iFreq) + Ha/rhoInterp(iFreq);
    QTilde = Qp/K_eq(iFreq) + Qa/KInterp(iFreq);
    CTilde = Cp*gamma_tilde(iFreq) + Ca*gamma_tildeInterp(iFreq);
    StiffTilde = [KTilde - omega(iFreq)^2*MTilde    -CTilde ;
        - CTilde.'          HTilde/omega(iFreq)^2 - QTilde];
    % Applying Dirichlet-type BCs
    sol = zeros(3*nTotNodes,1) ; % 3 Dofs per node
    sol(presBCDofs) = pres ; % enforcing pressures on left boundary
    sol(pres0BCDofs) = pres0 ; % enforcing pressures on top boundary
    F = F -StiffTilde*sol ;
    sol(freeDofs) = full(sparse(StiffTilde(freeDofs,freeDofs)) \ sparse(F(freeDofs))) ; % solving only for unconstrained cases
    %     % to get normal fluid displacements
    
    %% PRESSURE FIELD PLOTS
    %         PressureFieldToggle=0;
    if iFreq==10 && PressureFieldToggle==1
        keyboard
        % surface plot
        figure
        for theta=0%:(2*pi/60):4*pi
            clf
            
            hold on
            for iElem=1:nTotElem
                nodes=mesh.connect{iElem};
                x=mesh.coords(nodes,1);
                y=mesh.coords(nodes,2);
                z=real(sol(nodes)*exp(-1i*theta));
                
                patch(x,y,z,'FaceAlpha',.5)
                if mesh.matType(iElem)==2
                    patch('Faces',[1 2 3 4],'Vertices',mesh.coords(nodes,:),'EdgeColor','black','FaceColor','none','LineWidth',2,'FaceColor','black','FaceAlpha',.5)
                end
                
                
                
            end
            title([' t= ' num2str(theta/(2*pi*frequencies(iFreq))) ' s'])
            pause(0.001)
        end
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
    
    Px1_dof=2*length(mesh.coords)+x1_node; % dof corresponding to pressure at microphone 1
    Px2_dof=2*length(mesh.coords)+x2_node; % dof corresponding to pressure at microphone 2
    P1=sol(Px1_dof); % complex fourier pressure amplitude at microphone 1
    P2=sol(Px2_dof); % complex fourier pressure amplitude at microphone 2
    k=2*pi*frequences(iFreq)/c; % wave number in air
%     Rn=(P1*exp(-1i*k*x2)-P2*exp(-1i*k*x1) ) / (P1*exp(1i*k*x2)-P2*exp(1i*k*x1)) ; % -B/A LKKK
    Rn=(P1*exp(-1i*k*x2)-P2*exp(-1i*k*x1) ) / (-P1*exp(1i*k*x2)+P2*exp(1i*k*x1)); % B/A Luc
    SAC_2mic(iFreq) = 1-(abs(Rn))^2;
    
   
    %% Impedance p/u  method for calculating absorption coefficient
    
    Fg = StiffTilde * [sol(1:2*nTotNodes)/omega(iFreq).^2 ; sol(2*nTotNodes+1:end) ] ;
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
    
    figure(3)
    set(gca,'FontSize',16)
    hold on
    grid on
    
    plot(frequences,SAC_imp,'b--x','Linewidth',2,'Displayname','Impedance method') ; % FEM soln
    plot(frequences,SAC_2mic,'r--o','Linewidth',2,'Displayname','Two microphone method') ; % FEM soln
    xlabel('Frequency (Hz)')
    ylabel('\alpha - absorption coefficient')
    title('Biot-Helmholtz model')
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
