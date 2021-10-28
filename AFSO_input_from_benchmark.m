%% Script that creates parameters for AFSO problem from benchmark

global parameters

if ~exist('prob','var')
    warning('set the benchmark first')
    disp('prob=benchmarks(1)');
end

probID=prob.problemID;
if probID(11)==1
    skeletonmodel='BiotHelmholtz';
elseif probID(11)==0
    skeletonmodel='EquivalentFluid';
else
    error('Incorrect skeleton model')
end
frequencies=probID(5):probID(6):probID(7);

%% Inputs
plotToggle=0;
MeshplotToggle=0;
FEorder=1;
DisplayToggle=0;
pressurefieldToggle=0;
%% MEMOIZATION
% clearAllMemoizedCaches
% objfunc=memoize(@griddecoder1);
% objfunc.CacheSize=30;
% stat=stats(objfunc);
%% NO MEMOIZATION
% objfunc=(@griddecoder1);
%% Meshing Inputs
% common mesh


mesh.D=probID(3);
mesh.d=probID(4);
if ~isfield(mesh,'L')
mesh.L=0.05;
else
    warning('using existing L')
end
if ~isfield(mesh,'NELEMxL')
mesh.NELEMxL=10;
else
    warning('using existing NELEMxL')
end
mesh.NELEMxD=probID(1);
mesh.NELEMy=probID(2); % must be even for comparing full FE with symmetric FE
mesh.VF=probID(10);


if probID(8)==0
    mesh.bc='symmetric'; % 'full' assumes the mesh to be the whole impedance tube, 'symmetric' considers mesh to be only bottom half
elseif probID(8)==1
    mesh.bc='full';
end

if probID(9)==0
    mesh.sidebc='symmetric';
elseif probID(9)==1
    mesh.sidebc='fixed';
else
    error('Incorrect sidewall boundary condition')
end



%% MATERIAL CARD


%% meshing

[mesh.coords,mesh.connect,mesh.matType,mesh.type] = quad_mesher(mesh.L,mesh.D,mesh.d,mesh.NELEMxL,mesh.NELEMxD,mesh.NELEMy,mesh.bc);
mesh.domain=find(mesh.matType==2); 
switch mesh.bc
    case 'symmetric'
        mesh.N=mesh.NELEMxD*mesh.NELEMy; % for griddecoder1
        genomelength=mesh.N;
    case 'full'
        mesh.N=mesh.NELEMxD*mesh.NELEMy; % for griddecoder
        genomelength=mesh.N/2;
end


%% higher order mesh components
meshType = mesh.type;%'TRI' ;
if(strcmp(meshType,'QUAD'))
    Method = 'FEM2DQ4' ;
elseif(strcmp(meshType,'TRI'))
    Method = 'FEM2DT3' ;
end
if FEorder>1
    mesh = highOrderMesh(mesh,FEorder,Method) ; % For higher order meshing
    mesh.higherOrder=1;
end
mesh.connect = num2cell(mesh.connect,2); % problematic 


%% setting global parameters 
parameters{1}=mesh; % all the mesh and boundary condition information
parameters{2}=frequencies; % frequency vector
parameters{3}=plotToggle; % plots absorption curve 
parameters{4}=MeshplotToggle; % plots mesh
parameters{5}=FEorder; % order of finite element mesh
parameters{6}=DisplayToggle; % displays additional info
parameters{7}=prob.mat; % material properties
parameters{8}=skeletonmodel; % BiotHelmholtz or EquivalentFluid
parameters{9}=pressurefieldToggle; % plots pressure field after solving


%% sample genome set to fill the design domain fully with porous materials
genome=1*ones(genomelength,1);
mesh.matType(mesh.domain)=genome+1;

