% Initial DNA library used in first transcription
div=DivTrack(1000,1600,1e-10,'A',[0.2,0.8]);

fprintf('\nR1. Gel selection of cleaved\n');
div.measured('DNA for round 1 - BT382-A  E141',1000,1600);
div.T7(1300);
div.Select(true,.33);
div.measured('Post-gel',25,5000);
div.measured('Post-selection',240,133);
div.PCR(240,3400);

fprintf('\nR2. Gel selection of cleaved, no theo E148-E150\n');
div.measured('DNA for round 2 - BT383-B',1000,170);
div.T7(1720);
div.Select(true,0.25);
div.measured('pre-PCR DNA',1000,15);
div.PCR(239,3560);

fprintf('\nR3. +theo, no-cleavage, E153\n');
div.measured('DNA for round 3- BT392-A',400,100);
div.T7(1440);
div.measured('Post-RT',240,133);   % Where do these numbers come from?
div.Select(false,0.5);   % Guessing
div.PCR(239,2680);   % E154 - BT394

fprintf('\nR4 +theo, no-cleavage E155\n');
div.measured('DNA for round 4 - BT394-A',400,100);
%div.T7(30485*38/400);
div.T7(1158*1e-12/(div.volume*1e-6)*1e9);
div.measured('Post-RT',170,734);   % Where do these numbers come from? (were in E156 spreadsheet)
div.Select(false,0.25);   % Guessing?
div.PCR(240,2680);  % E156 - coincidence that same as R3 result
BT395=copy(div);

fprintf('\nR5a +theo, no-cleavage E157\n');
div.measured('DNA for round 5a - BT395-A',200,200);
div.T7(800*1e-12/(div.volume*1e-6)*1e9);
div.measured('Post-RT',170,1200);   % Where do these numbers come from? (were in E158 spreadsheet)
div.Select(false,0.10);   % Guessing?
div.PCR(240,2570);  % E156 - coincidence that same as R3 result
BT396=copy(div);

div=copy(BT395);
fprintf('\nR5b -theo, +cleavage E162\n');
div.measured('DNA for round 5b - BT395-A',200,200);
div.T7(278*1e-12/200e-6*1e9);
div.measured('RT-Input',100,2000);   
div.measured('Ligation product',400,84.6);
div.Select(true,0.041);   % Based on TRP
div.PCR(200,2880);
BT397=copy(div);

fprintf('\nR6b +theo no-cleavage E164\n');
div.measured('DNA for round 6b - BT397-B',200,200);
div.T7(277*1e-12/200e-6*1e9);
div.measured('Post-RT',100,1528);
div.Select(false,0.45);
div.PCR(200,2350);
BT398=copy(div);

fprintf('\nR7 -theo +cleavage IN PROCESS\n');
div.measured('DNA for round 7 - BT398-B',200,200);
div.T7(281.8*1e-12/200e-6*1e9);
div.measured('RT-Input',100,2000);   
%div.measured('Ligation product',?,?);
div.Select(true,0.5);   % Based on TRP
div.PCR(200,2600);
BT399=copy(div);

fprintf('\nR8 +theo -cleavage\n');
div.measured('DNA for round 8 - BT399-A',20,50);
div.T7(500);
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,15,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.Select(false,0.4);
div.measured('Part for PCR',25,div.conc);
div.PCR(100,450);

fprintf('\nR9 -theo +cleavage\n');
div.measured('DNA for round 9',15,div.conc/10);
div.T7(500);
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,15,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.Select(true,0.17);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.measured('Part for Ligation',69.0/3,div.conc);
div.measured('Ligation',div.volume*3,div.conc/3);
div.measured('Ligation Efficiency',div.volume, div.conc/2);
div.measured('Part for PCR',50,div.conc);
div.PCR(200,450);

fprintf('\nR10 +theo -cleavage\n');
div.measured('DNA for round 10 - BT401-B',15,50);
div.T7(500);
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,10,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.Select(false,0.4);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,450);

fprintf('\nR11 -theo +cleavage TODO\n');
div.measured('DNA for round 11',15,div.conc/12);
div.T7(500);
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,6,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.Select(true,0.17);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.measured('Part for Ligation',40.0/3,div.conc);
div.measured('Ligation',div.volume*3,div.conc/3);
div.measured('Ligation Efficiency',div.volume, div.conc/2);
div.measured('Part for PCR',25,div.conc);
div.PCR(100,450);

fprintf('\nR12 +theo -cleavage\n');
div.measured('DNA for round 12 - BT403-A',18,50);
div.T7(250);
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,11,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.Select(false,0.34);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,557);

fprintf('\nR13 -theo +cleavage\n');
div.measured('DNA for round 13',17,div.conc/12);
div.T7(447);
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,11,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.Select(true,0.40);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.measured('Part for Ligation',30.0/3,div.conc);
div.measured('Ligation',div.volume*3,div.conc/3);
div.measured('Ligation Efficiency',div.volume, div.conc/2);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,568);

fprintf('\n+theo -cleavage\n');
div.measured('DNA for round 14',18,50);
div.T7(250);   %  Unmeasured
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,11,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.Select(false,0.6);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,526);

fprintf('\n-theo +cleavage\n');
div.measured('DNA for round 15',17,div.conc/12);
div.T7(447);   % Unmeasured
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,11,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.Select(true,0.69);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.measured('Part for Ligation',30.0/3,div.conc);
div.measured('Ligation',div.volume*3,div.conc/3);
div.measured('Ligation Efficiency',div.volume, div.conc/2);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,529);

fprintf('\n+theo -cleavage\n');
div.measured('DNA for round 16',18,50);
div.T7(250);   %  Unmeasured
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,11,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.Select(false,0.42);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,444);

fprintf('\n-theo +cleavage\n');
div.measured('DNA for round 17',17,div.conc/12);
div.T7(447);   % Unmeasured
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,11,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.Select(true,0.48);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.measured('Part for Ligation',30.0/3,div.conc);
div.measured('Ligation',div.volume*3,div.conc/3);
div.measured('Ligation Efficiency',div.volume, div.conc/2);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,509);

fprintf('\n+theo -cleavage\n');
div.measured('DNA for round 18',18,50);
div.T7(250);   %  Unmeasured
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,11,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.Select(false,0.71);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,500);  % Unmeasuread

fprintf('\n-theo +cleavage\n');
div.measured('DNA for round 19',17,div.conc/12);
div.T7(447);   % Unmeasured
div.measured('Stop',div.volume*2,div.conc/2);
div.measured('Subsample RNA' ,11,div.conc);
div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
div.Select(true,0.69);
div.measured('Post-RT Dilution',div.volume*2,div.conc/2);
div.measured('Part for Ligation',30.0/3,div.conc);
div.measured('Ligation',div.volume*3,div.conc/3);
div.measured('Ligation Efficiency',div.volume, div.conc/2);
div.measured('Part for PCR',12.5,div.conc);
div.PCR(50,500);  % Unmeasuread

div.measured('Post Qiaquick+dilution',50,80);
for r=20:2:27
  fprintf('\n+theo -cleavage\n');
  div.measured(sprintf('DNA for round %d',r),10,div.conc*3.0/10);
  div.T7(75);   %  Unmeasured
  div.measured('Stop',div.volume*2,div.conc/2);
  div.measured('Subsample RNA' ,5,div.conc);
  div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
  div.measured('Post-RT Dilution',div.volume*3,div.conc/3);
  div.Select(false,0.5);  % Unmeasured
  div.measured('Part for PCR',25.0/4,div.conc);
  div.PCR(25,250); % Unmeasured
  div.measured('PCR Dilution',150,div.conc/6);
  
  fprintf('\n-theo +cleavage\n');
  div.measured(sprintf('DNA for round %d',r+1),10,div.conc*3.0/10);
  div.T7(250);   % Unmeasured
  div.measured('Stop',div.volume*2,div.conc/2);
  div.measured('Subsample RNA' ,5,div.conc);
  div.measured('Post-RT',div.volume*2,div.conc/2*0.8);
  div.Select(true,0.5);  % Unmeasured
  div.measured('Post-RT Dilution',div.volume*3,div.conc/3);
  div.measured('Part for Ligation',19.0/3,div.conc);
  div.measured('Ligation',div.volume*3,div.conc/3);
  div.measured('Ligation Efficiency',div.volume, div.conc/2);
  div.measured('Part for PCR',25.0/4,div.conc);
  div.PCR(50,250);  % Unmeasured
  div.measured('PCR Dilution',150,div.conc/6);
end
