% Track 1/28/15 large scale test of N5N60 TRP
% https://www.evernote.com/shard/s240/nl/4571194/6ee8e7ed-c546-4b93-bf82-87278a5219f5/
divsetup

goodFrac=1e-13*.171;    % Assume good molecules (i.e. aptamers) occur at the rate of 10^-13 of CLEAVABLE molecules
template=1-ceultramer^161;
templatep=(1-pagepurity)*template;

fprintf('\nPilot Run 1/28/15:\n');
div=DivTrack(1000,50,goodFrac,[bt88,templatep],'W',[0.2,0.8]);
div.T7(div.conc()*39.15); 
div.measured('Post-T7',49,39150);
div.Select(true,0.171);
div.measured('Post-PAGE',32,4231);
div.RT(1.0,bt600);
div.measured('Post-Ligation',704,30);

% Now run planned setup
fprintf('\nPlanned:\n');
div=DivTrack(1000,150,goodFrac,[bt88,templatep],'W',[0.2,0.8]);
div.T7(1958);
div.volume=50;
div.Select(true,0.171);
div.volume=50;
div.randchoose('Post-PAGE',.454);
div.RT(1.0,bt600);
div.volume=div.volume*12;
div.randchoose('Post-Ligation',.156*2);
% PCR using page-purified BT575
cycles=3;
pcrgain=2^cycles;
div.dilute(250/pcrgain,'Pre-PCR dilution');
div.PCR(cycles,[bt28,bt575p]);

div.randchoose('Use part for next round',2/pcrgain);
div.T7(1958);
