DivTrack - track diversity of a series of rounds of selection

Brent Townshend - 2/8/15

The DivTrack class tracks the number of molecules and diversity maintained over rounds of a selection experiment.  To do this, it makes several assumptions:

- the pool of molecules fall in 4 classes:
  - "good" molecules which switch;  cleavage is higher without target than with target, with all "good" molecules assumed to have the same clevage (e.g. 80% without, 20% with)
  - "bad" molecules whose cleavage rate is independent of the presence of target
  - "ragged" molecules, which have their 5' ends truncated.  There is a separate count maintained for ragged on the "left" vs "right" ends.  "left" is the sequence at the 5' end of the RNA, "right" is the 3' end of the RNA.   For example, the right end of the cDNA is actually the 5' end.   Ragged molecules are assumed to not cleave, and are handled differently for PCR amplification and transcription.  

- the pool of molecules initially contains a given fraction of molecules which are "good".  Each of the "good" molecules is assumed to initially have a unique sequence, but will later have multiple copies. The distribution of the number of occurences of a sequence is assumed to follow a truncated Poisson distribution (truncated in the sense that sequences occuring in 0 molecules are simply not in the mixture any more, so P(0)=0).  Thus, the probability of having k copies of a sequence is
      P(k)=q^k*exp(-q)/[(1-exp(-q))*k!]
      where q is the lambda of the Poisson distribution

- ragged molecules occur only due to the coupling efficiency (ce) during oligo synthesis.  The fraction of the molecules synthesized that are full length is ce^length. 

Several step of processing cause the pool of molecules to be sampled, with or without replacement.  For example, gel extraction results in some fraction of the molecules to be recovered, without replacement.  Transcription, however, samples with replacement -- the same molecule may be transcribed multiple times. The sampling steps result in changes to the distribution of "good" sequences. This distribution is assumed to be characterized by a single parameter, "kgood", which is the mean number of copies of each good sequence that is present.  Note that kgood will always be >=1.0 since only sequences which are present are considered.

Internally, the class keeps track of the total number of molecules in each class and the 'kgood'.  In addition, it uses a 'volume' property to be able to track the current volume and concentration of the pool.

The processing steps are implemented as the following core methods in the DivTrack class:

DivTrack(initvol, initconc, initfracgood,initfracragged,prefix,tgtCleave)
% initialize pool with given volume (in ul) and concentration (in nM)
% initfracgood is the fraction of the total pool that is "good"
% initfracragged is a two-element vector providing the fraction of the pool that is ragged on the left and right ends
% prefix is the prefix of the library;  'W','AW','BW','-'(cleaved);  just for keeping track
% tgtCleave is the cleavage without and with target for the "good" sequences

resample(obj,note,gain)
% Resample (with replacement) pool with given gain
% All of the molecules are sampled with the same gain, which can be any positive number
% note is a note to tie to this step

randchoose(obj,note,goodgain,badgain,raggedgain)
% Choose from the pool without replacement
% A different gain can be applied to each molecule class, each of which should be between 0 and 1

The above methods are used by several "helper" methods to simplify use:

Select(obj,keepCleave, cleaveRate)
% Choose molecules from the pool that either cleave or not (depending on keepCleave)
% The cleavage rate of the bad sequences is 'cleaveRate'.  Ragged molecules never cleave and good molecules cleave with tgtCleave

PCR(obj, cycles, primersRaggedFrac, isds)
% PCR amplify the pool with the given number of cycles 
% Assumes perfect amplication and uniform copying of all non-ragged input molecules
% The resulting pool will inherit the ragged fraction of the primers + the original unamplified ragged ones
% isds can be set to true to indicate the input is double-stranded; otherwise, if omitted, assumes input is single-stranded cDNA

T7(obj,rnaconc)
% Transcribe the pool ending with the given RNA concentration
% Assume molecules with ragged 5 end won't transcribe, but ragged 3' end ones are inherited
% After this method, the pool will be referring to the RNA produced only, not including the template

measured(obj,note,volume,conc)
% Update the pool to correspond to a measurement of an actual experiment
% Assumes that intervening steps caused loss of sample so applies randchoose with equal gain to all classes

