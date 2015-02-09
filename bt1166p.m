% Track 1/28/15 large scale test of N5N60 TRP
% https://www.evernote.com/shard/s240/nl/4571194/6ee8e7ed-c546-4b93-bf82-87278a5219f5/
ce=0.99;		% Coupling efficiency of oligo production (normal)
ceultramer=0.99; 	% Coupling efficiency of oligo production (Ultramers) -- adjusted to match observed qPCR-able frac
pagepurity=0.90;	% Purity after PAGE purification (IDT says ">90%", assume that 90% of bad products removed)
goodFrac=1e-13*.171;    % Assume good molecules (i.e. aptamers) occur at the rate of 10^-13 of CLEAVABLE molecules

% Fraction ragged for each oligo
% In most cases, if the oligo is shortened to less than 10 nt, it probably won't anneal anyway, 
%   so it doesn't carry forward the raggedness
bt600=1-ce^(28-10)/ce^5;	  % Will likely prime if it is at least 5 long, can be 10 short since PCR will be able to repair
bt28=1-ce^(21-10);	% Assume 10 short will still get PCR'ed by B/X
bt88=1-ce^(23-10);	% Assume 10 short will still get transcribed
bt575=1-ce^29;		% Any shortening will affect subsequent transcription
bt575p=(1-pagepurity)*bt575;	
bt1076=1-ce^62;		  % Extension oligos need to be full-length to work at all
bt1076p=(1-pagepurity)*bt1076;
bt654=1-ce^(43-16)/ce^8;   % Will likely splint extension if hasn't been shortened by more than 16 (at least 1 overhang),
bt654p=(1-pagepurity)*bt654;  % Most will probably still splint
			  %   and will likely block cDNA if it is 8 or more long (full=43)
template=1-ceultramer^161;
templatep=(1-pagepurity)*template;

fprintf('Fraction of product from each oligo that will be ragged on 5'' end:\n');
fprintf('Template=%.2f/%.2f, BT600=%.2f, BT28=%.2f, BT88=%.2f, BT575=%.2f/%.2f, BT1076=%.2f/%.2f, BT654=%.2f/%.2f\n', template, templatep, bt600,bt28,bt88,bt575,bt575p,bt1076,bt1076p,bt654,bt654p);

fprintf('\nPilot Run 1/28/15:\n');
div=DivTrack(1000,50,goodFrac,[bt88,templatep],'W',[0.2,0.8]);
div.T7(div.conc()*39.15); 
div.measured('Post-T7',49,39150);
div.Select(true,0.171);
div.measured('Post-PAGE',32,4231);
div.RT(1.0,bt600);
div.measured('Post-Ligation',704,30);
% 3-cycle PCR using page-purified BT575
div.PCR(div.volume,div.conc*2^3,[bt28,bt575p]);
setfig('bt1166p'); clf;
div.plothistory(); 
suptitle(sprintf('Observed BT1166p $%.0f',div.cumcost));

% Now run planned setup
fprintf('\nPlanned:\n');
div=DivTrack(1000,150,goodFrac,[bt88,templatep],'W',[0.2,0.8]);
div.T7(1958);
div.volume=49;
div.Select(true,0.171);
div.volume=32;
div.randchoose('Post-PAGE',.454);
div.RT(1.0,bt600);
div.volume=704;
div.randchoose('Post-Ligation',.156*2);
% 3-cycle PCR using page-purified BT575
div.PCR(div.volume,div.conc*2^3,[bt28,bt575p]);
setfig('bt1166p'); clf;
div.plothistory(); 
suptitle(sprintf('Observed BT1166p $%.0f',div.cumcost));
