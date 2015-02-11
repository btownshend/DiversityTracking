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
div.plotdist('Trial');

% Now run planned setup
fprintf('\nPlanned:\n');
div=DivTrack(1200,150,goodFrac,[bt88,templatep],'W',[0.2,0.8]);
div.T7(1958);
div.volume=50;
div.Select(true,0.171);
div.volume=50;
div.randchoose('Post-PAGE',.454);
div.dilute(0.8*3150);	% Dilute to 80% of Omniscript max capacity
div.RT(1.0,bt600);
div.dilute(min(div.conc/12,100));	% Keep below splint+oligo conc, but with at least 12x dilution
div.randchoose('Post-Ligation',.24*1.84);	% Assuming input in pilot run was limited to 125nM extension concentration, efficiency was 30/125=.24
                                                % Using PAGE oligos improves this by 1.84x
div.plotdist('R1 PostLig');

% PCR using page-purified BT575
cycles=4;
pcrgain=2^(cycles-1);
div.dilute(1000/pcrgain,'Pre-PCR dilution');
div.PCR(cycles,[bt28,bt575p]);

div.volume=div.volume*3/div.kgood;
div.randchoose('Use part for next round',4/div.kgood);  % Keep 2 copies/sequence
div.dilute(200);	% Tune to not lose any more diversity
div.T7(1958);
div.volume=50;
div.Select(true,0.3);   % Assume 30% cleavage
div.volume=50;
div.randchoose('Post-PAGE',.454);
div.dilute(0.8*3150);
div.RT(1.0,bt600);
div.dilute(min(div.conc/12,100));	% Keep below splint+oligo conc, but with at least 12x dilution
div.randchoose('Post-Ligation',.24*1.84);
div.plotdist('After R2 Ligation');
% PCR using page-purified BT575
cycles=3;
pcrgain=2^(cycles-1);
div.dilute(1000/pcrgain,'Pre-PCR dilution');
div.PCR(cycles,[bt28,bt575p]);
div.plotdist('After R2 PCR');
