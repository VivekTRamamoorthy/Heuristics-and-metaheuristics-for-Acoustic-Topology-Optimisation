%% benchmark problem instances
function out=benchmarks(problemno)
%%
bmnames={'LKKKcoarse','building','low-f-high-$\sigma$','Automotive','highfreq','LKKKfine','500Hzproblem'};
benchmarkproblems=...
[10 10 0.135 0.054  100   100   1500      0 0 1.0  1 1; % LKKK 1
15  10 0.045 0.100  100   100   1500      0 0 1.0  1 2; % Building 2
50  20 0.100 0.100  50    50    500       0 0 1.0  1 3; % low frequency high resistivity 3
10  10 0.020 0.100  100   100   1500      0 0 1.0  1 2; % auto 4
10  10 0.020 0.100  2000  1000  5000      0 0 1.0  1 2; % high frequency application 5
50  20 0.135 0.054  100   100   1500      0 0 1.0  1 2; % LKKK fine mesh melamine  6
10  5  0.135 0.054  500   500   500       0 0 1.0  1 2; % fast problem 500 hz 7
50  20 0.135 0.054  100   100   1500      0 1 1.0  1 2; % LKK actual problem 8
10  5  0.135 0.054  500   500   500       0 0 1.0  0 2; % Prob inst 7 but EF at just 500 9
50  20 0.135 0.054  500   500   500       0 0 1.0  1 2; % Prob for testing pressure field algorithm 10 BH
50  20 0.135 0.054  500   500   500       0 0 1.0  0 2; % Prob for testing pressure field algorithm 11 EF
50  20 0.34 0.054  500   500   500       0 0 1.0  0 2;]; % Prob for testing pressure field algorithm 12 EF
%nelx nely D  y     f1    fstep f2     full fixedsidewall volfrac %BH %mat


% 
% 10  10 0.135 0.054 100  1    100      0 0 0.5 1;
% 10  10 0.135 0.054 200  1    200      0 0 0.5 1;
% 10  10 0.135 0.054 500  1    500      0 0 0.5 1;
% 10  10 0.135 0.054 1000 1    1000     0 0 0.5 1;
% 10  10 0.135 0.054 2000 1    2000     0 0 0.5 1;
% 10  10 0.135 0.054 5000 1    5000     0 0 0.5 1;
% 
% 50  20 0.135 0.054 100  100  1500     0 0 0.5 1];

%%
% LKKKprobID=[ 50 20 0.135 0.054 100 100 1500  0 1 0.5 1 1];
% standardprobID=[10 10 0.135 0.054 50 500 19  0 0 0.5 1 2];
%%


problemID=benchmarkproblems(problemno,:);
% problemID=[ (mesh.NELEMxD) (mesh.NELEMy) (mesh.D) (mesh.d) frequencies(1) frequencies(end) length(frequencies)  strcmp(mesh.bc,'full') strcmp(mesh.sidebc,'fixed')  ];
% num2str(problemID')
% probleminstance=[ num2str(mesh.NELEMxD) ' x ' num2str(mesh.NELEMy) ' grid ' num2str(mesh.D) ' m x' num2str(mesh.d) ' m ' num2str(frequencies(1)) ' to ' num2str(frequencies(end)) ' Hz '  num2str(length(frequencies)) ' steps ' mesh.bc ' problem ' mesh.sidebc ' sidewall' ];

% probleminstance=[ num2str(problemID(1)) ' x '  num2str(problemID(2)) ' grid, from ' ...
%     num2str(problemID(3)) ' m x' num2str(problemID(4)) ' m ' ...
%     num2str(problemID(5)) ' in steps of ' num2str(problemID(6)) ' upto '  num2str(problemID(7)) ' Hz '  num2str(problemID(8)) ...
%     ' fullproblem '  num2str(problemID(9)) ' fixedsidewall ' ...
%     num2str(problemID(10)) ' volfrac ' num2str(problemID(11)) ' poroelastic'];

probleminstance={['Grid: ' num2str(problemID(1)) ' x '  num2str(problemID(2)) ];
    ['Dimensions: '  num2str(problemID(3)) ' m x' num2str(problemID(4)) ' m ' ];
    ['Frequencies : ' num2str(problemID(5)) ' : ' num2str(problemID(6)) ' : '  num2str(problemID(7)) ' Hz '];
    ['Full problem? (1-y,0-n): ' num2str(problemID(8)) ];
    ['Sidewall fixed? (1-y,0-n): '  num2str(problemID(9)) ]; 
    ['Vol frac : ' num2str(problemID(10))]; 
    ['BH poroelastic? : (1-y,0-n)' num2str(problemID(11)) ];
    ['Material no: ' num2str(problemID(12))]};

out.problemID=problemID;
% disp('Problem Instance:')
% disp(probleminstance)
out.problem=probleminstance;

materialno=problemID(12);
%% Material selection
mat=matdatabase(materialno);

% disp('Material properties')
% disp(materialcard)
out.mat=mat;
%% Design of experiments
% standardprobID=[ 10 10 0.135 0.054 50 25 500  0 0 0.5 1];
% meshes=[10 10;
%     20 10;
%     50 20];
% dimensions=[0.06 0.054;
%     0.135 0.054;
%     0.27 0.054];
% frequencies= [100 100 1500;
%     50 25 500 ;
%     50 25 250;
%     300 100 1500;
%     2000 1000 5000;
%     100 1 100
%     500 1 500
%     1000 1 1000
%     2000 1 2000];
% fullprob=[0;1];
% sidebcfixed=[0;1];
% volfracs=[0.5;0.25;0.1];
% skeletonmodel=[0;1];
% 
% LKKKprobID=[meshes(3,:) dimensions(2,:) frequencies(1,:) 0 1 0.5 1];



%% printing benchmarks in latex table format
tabtoggle=0;
if tabtoggle
    %%
    table={['\begin{table}'];
        ['\centering'];
        [' \caption{List of Benchmark problem instances} '];
        ['   \begin{tabular}{c c c c c c c c c c c}'];
        ['   \hline  NELEMx & NELEMy & Length & Height & $f_{min}$ &$f_{step}$ & $f_{max}$ & fullFE & BC fixed & $V_f$ & $BH$ \\ \hline']};
    
    for i=1:length(benchmarkproblems(:,1))
        string=['    '];
        for j=1:length(benchmarkproblems(1,:))-1
            string=[string num2str(benchmarkproblems(i,j)) ' & '];
        end
        string=[string num2str(benchmarkproblems(i,j+1)) ' \\ '];
        table{length(table)+1}=[string ' '];
        
    end
    table(end+1:end+4)={['   \hline '];
        ['   \end{tabular} '];
        ['   \label{tab:benchmarks} '];
        ['\end{table}']};
    
    
    cells2text(table,'table.out')
    
end

%% printing benchmarks in latex table format
tabshorttoggle=0;



if tabshorttoggle
    %%
    table={['\begin{table}'];
        ['\centering'];
        %[' \caption{List of Benchmark problem instances} '];
        ['   \begin{tabular}{|c c c c c c c c c c|}']};
       % ['   \hline  NELEMx & NELEMy & Length & Height & $f_{min}$ &$f_{step}$ & $f_{max}$ & fullFE & BC fixed & $V_f$ & $BH$ \\ \hline']};
    
    
    
    colnames={ ' NELEMx ' ' NELEMy ' 'Length ' ' Height ' ' $f_{min}$  '  ' $f_{step}$   ' ' $f_{max}$  ' ' fullFE  ' 'BC fixed  ' ' $V_f$  ' ' $BH$ ' ' Mat no'};
    colstoinclude=[1 2 3 4 5 6 7 10 12]; 
    
        % col names generation
        string=['   \hline Name &'];
        for j=1:length(colstoinclude)-1
            string=[string colnames{colstoinclude(j)} ' & '];
        end
        string=[string colnames{colstoinclude(j+1)} ' \\ '];
        table{length(table)+1}=[string ' \hline'];
        
        % cols generation

        for i=1:length(benchmarkproblems(:,1))
            string=['    ' bmnames{i} ' & '];
            for j=1:length(colstoinclude)-1
                string=[string num2str(benchmarkproblems(i,colstoinclude(j))) ' & '];
            end
            string=[string num2str(benchmarkproblems(i,colstoinclude(j+1))) ' \\ '];
            table{length(table)+1}=[string ' '];
            
        end
        table(end+1:end+4)={['   \hline '];
            ['   \end{tabular} '];
            ['   \label{tab:benchmarks} '];
            
            ['\end{table}']};
        
        
        cells2text(table,'tableshort.out')
        open tableshort.out
        
end
