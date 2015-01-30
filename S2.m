% Initial DNA library used in first transcription
geleff=0.5;
coleff=0.8;
ligeff=0.25;
fullfrac=1;
pcrgain=8;

fprintf('\nR1a. Gel selection of cleaved from PAGE-purified input\n');
div=DivTrack(1000,50,1e-10,'W',[0.2,0.8]);
div.T7(1500);
div.diluteToVolume(1012);
div.Select(true,.3);
div.randchoose('gel efficiency',geleff);
div.diluteToVolume(100,'Zymo-Clean');
div.diluteToVolume(div.volume/0.7,'RT');
div.randchoose('RT efficiency',0.8);
div.diluteToVolume(1.1*div.volume,'Ligation');
div.randchoose('lig efficiency',ligeff);
div.dilute(250/pcrgain,'PCR Input dilution');
div.PCR(div.volume,div.conc*pcrgain);

return 

fprintf('\nR1b. Ligation selection of cleaved from PAGE-purified input\n');
div=DivTrack(500,200,1e-10,'A',[0.2,0.8]);
div.T7(1510);
div.randchoose('lig efficiency',ligeff,ligeff);
div.Select(true,.3);
div.dilute(250/2^4);
div.PCR(div.volume,250);

fullfrac=0.15;
fprintf('\nR1c. Gel selection of cleaved from non-purified input with %.2f full-length\n',fullfrac);
div2=DivTrack(500,200,1e-10*fullfrac,'A',[0.2,0.8]);
gain=16;  % Amount of gain of full-length
div2.dilute(250/gain);
adj=gain/(fullfrac*gain+(1-fullfrac));
fullfrac=fullfrac*adj;
div2.fracgood=div2.fracgood*adj;
div2.randchoose('PCR of full-length',fullfrac+(1-fullfrac)/gain);
div2.PCR(div2.volume,250);
div2.dilute(1510);
div2.T7(1510);
div2.Select(true,.3*fullfrac);
div2.randchoose('gel efficiency',geleff,geleff);
div2.randchoose('lig efficiency',ligeff,ligeff);
div2.dilute(250/2^4);
div2.PCR(div2.volume,250);
