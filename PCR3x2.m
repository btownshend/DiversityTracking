% BT1175 - 12%
divsetup
fracGood=.12;
badCleavage=0.171;

div=DivTrack(480,125/fracGood,1e-13*badCleavage,[0,1-fracGood],'W',[0.2,0.8]);
div.PCR(div.volume,div.conc*2,[bt88,bt575]);
div.resample('keep 1/2',0.5);
div.PCR(div.volume,div.conc*2,[bt88,bt575]);
div.resample('keep 1/2',0.5);
div.PCR(div.volume,div.conc*2,[bt88,bt575]);
