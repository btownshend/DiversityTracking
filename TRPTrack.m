% TRPTrack - track diversity of a series of rounds of TRP

classdef TRPTrack < matlab.mixin.Copyable
  properties
    volume; % Volume in ul
    ngood;  % Number of good switch molecules in mixture
    nbad;   % Number of other molecules in mixture (i.e. non-switches)
    initfracgood;
    initngood;
    tgtCleave;
    prefix;
    history;
    cumcost;
    randtargets;   % Number of unique random numbers to use for good seqs
    goodseqs;  % Each value in this vector represents a unique sequence, keeps empirical distribution
    poolCleavage;  % Empirical distribution of pool cleavages
    poolFrac;	   % Corresponding prob of each pool cleavage
  end

  methods(Static)
    function m=moles(vol,conc)
      m=vol*1e-6*conc*1e-9*6.022e23;
    end
  end
  
  
  methods
    function obj=TRPTrack(initvol, initconc, initfracgood,prefix,tgtCleave,poolCleavage)
    % initialize pool with given volume (in ul) and concentration (in nM)
    % initfracgood is the fraction of the total pool that is "good"
    % prefix is the prefix of the library;  'W','AW','BW','-'(cleaved);  just for keeping track
    % tgtCleave is the cleavage with and without target for the "good" sequences      
    % poolFitness is an empirical distribution of the cleavage of the pool, each value is equally probable
      obj.volume=initvol;
      obj.ngood=obj.moles(initvol,initconc)*initfracgood;
      obj.nbad=obj.moles(initvol,initconc)*(1-initfracgood)
      obj.initfracgood=initfracgood;
      obj.initngood=obj.ngood;
      obj.randtargets=1000;
      obj.goodseqs=1:obj.randtargets;
      obj.prefix=prefix;
      obj.tgtCleave=tgtCleave;
      obj.history=[];
      obj.setPoolCleavage(poolCleavage);
      obj.poolFrac=ones(size(poolCleavage))/length(poolCleavage);
      fprintf('Initial prefix=%s, Target cleavage=[%.2f, %.2f], Initial good=%.3g bad=%.2g\n',obj.prefix,obj.tgtCleave,obj.ngood, obj.nbad);
      obj.printdiv('Initial');
      obj.cumcost=0;
    end

    function setPoolCleavage(obj,poolCleavage)
    % Update pool cleavage to given setting
      obj.poolCleavage=poolCleavage;
      obj.poolFrac=ones(size(poolCleavage))/length(poolCleavage);
    end

    function k=kgood(obj)
    % Compute the mean number of copies of each good sequence
      k=length(obj.goodseqs)/length(unique(obj.goodseqs));
    end
    
    function g=divtarget(obj)
    % Compute the fraction of the diversity target we are at
      g=length(unique(obj.goodseqs))/obj.randtargets * obj.initngood;
    end
    
    function t=total(obj)
      t=obj.ngood+obj.nbad;
    end

    function f=fracgood(obj)
      f=obj.ngood/obj.total();
    end
    
    function c=conc(obj)
    % Concentration in nM
      c=obj.total()/(obj.volume*1e-6)/6.022e23*1e9;
    end
        
    function printdiv(obj,note)
      fprintf('%-50.50s   %2s %5.0ful %6.1fnM Total=%7.2g Enrich=%7.2g kgood=%8.2f goodseqs=%5.3g cumcost=$%3.0f\n',note,obj.prefix, obj.volume,obj.conc(),obj.total(),obj.fracgood()/obj.initfracgood, obj.kgood(), obj.divtarget(), obj.cumcost);
      obj.history=[obj.history,struct('ngood',obj.ngood,'bad',obj.nbad,'divtarget',obj.divtarget,'note',note)];
    end

    function plothistory(obj)
      subplot(211);
      h=semilogy([obj.history.bad],'-o');
      hold on;
      set(gca,'XTick',1:length(obj.history));
      set(gca,'XTickLabel',{});
      legend(h,{'bad'});
      
      subplot(212);
      h=semilogy([obj.history.ngood],'-o');
      hold on;
      h(2)=semilogy([obj.history.divtarget],'-o');
      legend(h,{'good','good sequences'});
      set(gca,'XTick',1:length(obj.history));
      set(gca,'XTickLabel',{obj.history.note});
      set(gca,'XTickLabelRotation',15);
      c=axis;
      c(3)=0.1;
      axis(c);
    end
    
    function plotdist(obj,ti)
    % Plot distribution of counts of sequences
      if nargin<2
        ti='Divtrack.dist';
      end
      cnt=hist(obj.goodseqs,1:max(obj.goodseqs));
      [dist,n]=hist(cnt,0:max(cnt));
      dist=dist/sum(dist)*obj.divtarget();
      setfig(ti);clf;
      plot(n(n>0),dist(n>0),'o-');
      xlabel('Number of copies of sequence');
      ylabel('Number of sequences');
      title(sprintf('%s: Total of %.2f sequences with mean of %.2f copies',ti, obj.divtarget(), obj.kgood()));
    end
    
    function resample(obj,note,gain)
    % Resample (with replacement) pool with given gain
      obj.ngood=obj.ngood*gain;
      obj.nbad=obj.nbad*gain;
      obj.goodseqs=randsample(obj.goodseqs,round(gain*length(obj.goodseqs)),true);
      obj.printdiv(sprintf('%s(gain=%.3f)',note,gain));
    end

    function dilute(obj,tgtConc,note)
      obj.volume=obj.volume*obj.conc()/tgtConc;
      if nargin<3
        note='Dilute';
      end
      obj.printdiv(note);
    end
    
    function diluteToVolume(obj, volume, note)
      if nargin<3
        note='Dilute';
      end
      obj.volume=volume;
      obj.printdiv(note);
    end
    
    function randchoose(obj,note,goodgain,badgain)
    % Choose from the pool without replacement
    % A different gain can be applied to each molecule class, each of which should be between 0 and 1
      if nargin<4
        badgain=goodgain;
      end
      if badgain>1.0 || badgain<0.0 || goodgain>1.0 || goodgain<0.0
        error('Bad gains: %f, %f, %f\n', goodgain, badgain);
      end
      % Assume that there are exactly k copies of each sequence (TODO: does this reasonably approximate?)
      obj.ngood=obj.ngood*goodgain;
      obj.nbad=obj.nbad*badgain;
      obj.goodseqs=randsample(obj.goodseqs,round(goodgain*length(obj.goodseqs)),false);
      if goodgain==badgain
        obj.printdiv(sprintf('%s(gain=%.3f)',note,goodgain));
      else
        obj.printdiv(sprintf('%s(gain=%.3fG,%.3fB)',note,goodgain,badgain));
      end
    end

    function gamma=findGamma(obj,keepCleave,retentionFraction)
    % Find scale factor for non-cleavers that results in given retentionFraction
    % Assume (1-cleavage) is raised to some power 
      if keepCleave
        % i.e. retain(i)=1-(1-poolCleavage(i))^gamma;   sum(retain(i)*poolFrac(i))=retentionFraction
        gamma=fminsearch(@(z) (sum((1-(1-obj.poolCleavage).^z).*obj.poolFrac)-retentionFraction)^2,1.0);
      else
        % Find scale factor for non-cleavers that results in given retentionFraction
        % i.e. retain(i)=(1-poolCleavage(i))^gamma;   sum(retain(i)*poolFrac(i))=retentionFraction
        gamma=fminsearch(@(z) (sum((1-obj.poolCleavage).^z.*obj.poolFrac)-retentionFraction)^2,1.0);
      end
    end
      
    function Select(obj,keepCleave,gamma)
    % Choose molecules from the pool that either cleave or not (depending on keepCleave)
    % The cleavage rate of the bad molecules is based on the empirical distribution poolCleavage and poolFrac.  
    % Cleavage values are mapped to values that result in the observed retentionFraction
      if keepCleave
        retain=1-(1-obj.poolCleavage).^gamma;
        tgtRetain=1-(1-obj.tgtCleave(2))^gamma;
        meanRetain=sum(retain.*obj.poolFrac);
        obj.prefix='-';
        obj.randchoose(sprintf('Select(Clv,%.1f)',gamma),tgtRetain,meanRetain);
      else
        obj.prefix='W';
        fprintf('Gamma=%.1f\n',gamma);
        retain=(1-obj.poolCleavage).^gamma;
        tgtRetain=(1-obj.tgtCleave(1))^gamma;
        meanRetain=sum(retain.*obj.poolFrac);
        obj.randchoose(sprintf('Select(Unclv,%.1f)',gamma),tgtRetain,meanRetain);
      end
      obj.poolFrac=obj.poolFrac.*retain;
      obj.poolFrac=obj.poolFrac/sum(obj.poolFrac);
    end
    
    function T7(obj,rnaconc)
    % Transcribe the pool ending with the given RNA concentration
    % After this method, the pool will be referring to the RNA produced only, not including the template
      obj.cumcost=obj.cumcost+0.079*obj.volume;   % Price/ul of rx:  T7: .050 NTP .016, SuperaseIn .013
      gain=rnaconc*obj.total()/(obj.conc()*(obj.ngood+obj.nbad));
      obj.resample('T7',gain);
    end

    function BeadSep(obj, beadConc, efficiency)
    % Bead separation with given bead conc (mg/ml) and efficiency
      obj.randchoose('Bead separation',efficiency);
      cost=1648/(100e3);   % Cost of 1ug of beads  based on $1648/(10ml * 10mg/ml) 
      obj.cumcost=obj.cumcost+obj.volume*beadConc*cost;
    end
    
    function RT(obj, efficiency)
    % Reverse transcription
    % bad, good, all get same gain of efficiency
      obj.cumcost=obj.cumcost+.236*obj.volume;
      % Apply the gain
      obj.randchoose('RT',efficiency);
    end
    
    function ligate(obj, prefix, efficiency)
      obj.prefix=prefix;
      obj.randchoose(sprintf('ligate(%s)',prefix),efficiency);
      obj.cumcost=obj.cumcost+.0002*obj.volume;
    end
    
    function gain=PCR(obj,ncycles,isds)
    % PCR amplify the pool to the given final volume and concentration
    % Assumes perfect amplication and uniform copying of all input molecules
    % isds can be set to true to indicate the input is double-stranded; otherwise, if omitted, assumes input is single-stranded cDNA
      if nargin<4
        isds=false;
      end
      
      % Initial S1 to forward strand (same direction as RNA), S2 as reverse strand
      if isds
        S1=struct('nbad',obj.nbad,'ngood',obj.ngood,'goodseqs',obj.goodseqs);
      else
        S1=struct('nbad',0,'ngood',0,'goodseqs',[1]);  % Keep 1 sequence to avoid errors later
      end
      S2=struct('nbad',obj.nbad,'ngood',obj.ngood,'goodseqs',obj.goodseqs);

      for i=1:ncycles
        % Update
        %fprintf('S1=%.2g %.2g S2=%.2g %.2g \n', S1.nbad, S1.ngood, S2.nbad, S2.ngood);
        S1new=struct('nbad',S1.nbad+S2.nbad,...
                     'ngood',S1.ngood+S2.ngood,...
                     'goodseqs',[S1.goodseqs,randsample(S2.goodseqs,round(length(S2.goodseqs)))]);
        S2new=struct('nbad',S2.nbad+S1.nbad,...
                     'ngood',S2.ngood+S1.ngood,...
                     'goodseqs',[S2.goodseqs,randsample(S1.goodseqs,round(length(S1.goodseqs)))]);
        S1=S1new;
        S2=S2new;
      end
      obj.goodseqs=S2.goodseqs;   % Reverse-complement strand is template for future transcription
      fprintf('S1=%.2g %.2g  S2=%.2g %.2g\n', S1.nbad, S1.ngood, S2.nbad, S2.ngood);
      % No change to number of good sequences since sampling is uniform
      oldtotal=obj.total;
      % During last cycle, assume everything got extended so all parameters are aveages
      obj.ngood=mean([S1.ngood,S2.ngood]);
      obj.nbad=mean([S1.nbad,S2.nbad]);
      gain=obj.total/oldtotal;
      obj.cumcost=obj.cumcost+obj.volume*263/250/50;  % Kapa is $263/250U, uses 1U/50ul
      obj.printdiv(sprintf('PCR (gain=%.3f)',gain));
    end

    function measured(obj,note,volume,conc)
    % Update the pool to correspond to a measurement of an actual experiment
    % Assumes that intervening steps caused loss of sample so applies randchoose with equal gain to all classes
      obj.volume=volume;
      gain=obj.moles(volume,conc)/obj.moles(obj.volume,obj.conc());
      obj.randchoose(note,gain);
    end

    function keeppart(obj,note,volume)
    % Only keep given volume
      gain=volume/obj.volume;
      obj.volume=volume;
      obj.randchoose(note,gain);
    end
  end
end
