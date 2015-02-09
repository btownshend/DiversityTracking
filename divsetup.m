% Setup common parameters for diversity tracking

ce=0.99;		% Coupling efficiency of oligo production (normal)
ceultramer=0.99; 	% Coupling efficiency of oligo production (Ultramers) -- adjusted to match observed qPCR-able frac
pagepurity=0.90;	% Purity after PAGE purification (IDT says ">90%", assume that 90% of bad products removed)

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

fprintf('Fraction of product from each oligo that will be ragged on 5'' end:\n');
fprintf('Template=%.2f/%.2f, BT600=%.2f, BT28=%.2f, BT88=%.2f, BT575=%.2f/%.2f, BT1076=%.2f/%.2f, BT654=%.2f/%.2f\n', template, templatep, bt600,bt28,bt88,bt575,bt575p,bt1076,bt1076p,bt654,bt654p);

