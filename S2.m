% Initial DNA library used in first transcription
coleff=0.8;
ligeff=0.5;
geleff=0.464;
rteff=0.9;
fullfrac=1;
pcrgain=8;
ce=0.99;
lib=161;
badCleavage=0.171;
goodFrac=badCleavage/1e13;  % Want initial pool to have 10^-13 relative to CLEAVABLE sequences
pagepurity=0.92;
pageoligos=false;

for i=1:5
  fprintf('\n');
  if i>=4
    pageoligos=true;
  end
  % Compute fraction of results that will be ragged when using each oligo
  if pageoligos
    bt600=1-pagepurity;
    bt30=1-pagepurity;
    bt88=1-pagepurity;
    bt575=1-pagepurity;
    bt1076=1-pagepurity;
    bt654=(1-pagepurity)/2;  % Most will probably still splint
  else
    bt600=1-ce^((28-5)-10);
    bt30=1-ce^(45-10);
    bt88=1-ce^(23-10);
    bt575=1-ce^(29-10);
    bt1076=1-ce^62;
    bt654=1-ce^((43-16)-8);
  end
  fprintf('Fraction of product from each oligo that will be ragged on 5'' end:\n');
  fprintf('BT600=%.2f, BT30=%.2f, BT88=%.2f, BT575=%.2f, BT1076=%.2f, BT654=%.2f\n', bt600,bt30,bt88,bt575,bt1076,bt654);
  if i==1
    ti=sprintf('R1a. Gel selection of cleaved from PAGE-purified inputs\n');
    fprintf(ti);
    vol=1000; % High enough to keep kgood>=1.5
             % Start with at least 1.5 copies of good seqs
    div=DivTrack(vol,150,goodFrac,[1-pagepurity,1-pagepurity],'W',[0.2,0.8]);
  elseif i==2
    ti=sprintf('R1b. Gel selection of cleaved from non-PAGE-purified inputs\n');
    fprintf(ti);
    div=DivTrack(7000,400,goodFrac,[bt88,1-ce^lib],'W',[0.2,0.8]);
  elseif i==3
    ti=sprintf('R1c. Gel selection of cleaved from PCR products of non-PAGE-purified inputs\n');
    fprintf(ti);
    div=DivTrack(4000,105,goodFrac,[0.0,1-ce^lib],'W',[0.2,0.8]);
    div.PCR(div.volume,div.conc*8,[bt88,bt575]);
    div.volume=div.volume/4; div.randchoose('Part of PCR',0.25);
    div.diluteToVolume(2000,'T7 input dilution');
  elseif i==4
    ti=sprintf('R1d. Gel selection of cleaved from PCR products of non-PAGE-purified inputs using PAGE-oligos\n');
    fprintf(ti);
    div=DivTrack(2500,155,goodFrac,[0.0,1-ce^lib],'W',[0.2,0.8]);
    div.PCR(div.volume,div.conc*4,[bt88,bt575]);
    div.volume=div.volume/2; div.randchoose('Part of PCR',0.5);
    div.diluteToVolume(900,'T7 input dilution');
  elseif i==5
    ti=sprintf('R1e. Gel selection of cleaved from PCR products of PAGE-purified inputs using PAGE-oligos\n');
    fprintf(ti);
    div=DivTrack(1500,65,goodFrac,[0.0,1-pagepurity],'W',[0.2,0.8]);
    div.PCR(div.volume,div.conc*4,[bt88,bt575]);
    div.volume=div.volume/2; div.randchoose('Part of PCR',0.5);
    div.diluteToVolume(370,'T7 input dilution');
  end
  div.T7(1500);   % Assumed yield from T7 is 1.5uM, independent of template concentration
  
  div.Select(true,badCleavage);
  div.randchoose('gel efficiency',geleff);
  div.diluteToVolume(100,'Zymo-Clean concentrate');
  div.diluteToVolume(div.volume/0.7,'RT input dilution');
  div.RT(rteff,bt600);
  div.diluteToVolume(1.1*div.volume,'Ligation input dilution');
  div.ligate('BW',ligeff*(1-bt1076)*(1-bt654));
  div.dilute(250/pcrgain,'PCR Input dilution');
  div.PCR(div.volume,div.conc*pcrgain,[bt30,bt575]);
  setfig(sprintf('R1-%d',i)); clf;
  div.plothistory(); suptitle(sprintf('%s $%.0f/lib',ti,div.cumcost));
  fprintf('Final fraction good sequences = 1 in %.2g\n',div.total()/div.kgood);
  %fprintf('Total cost for 10 libraries:  $%.0f\n', div.cumcost*10);
end
