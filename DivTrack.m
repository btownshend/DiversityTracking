% DivTrack - track diversity of a series of rounds of TRP

classdef DivTrack < matlab.mixin.Copyable
  properties
    volume; % Volume in ul
    ngood;  % Number of good switch molecules in mixture
    nragged;	% Number with incorrect ends (due to short oligos) - first of the two elements refers to 5' end of RNA (or corresponding end of DNA)
    nbad;   % Number of other molecules in mixture (i.e. non-switches)
    initfracgood;
    initngood;
    tgtCleave;
    prefix;
    history;
    cumcost;
    goodseqs;  % Each value in this vector represents a unique sequence, keeps empirical distribution
  end

  methods(Static)
    function m=moles(vol,conc)
      m=vol*1e-6*conc*1e-9*6.022e23;
    end
  end
  
  
  methods
    function obj=DivTrack(initvol, initconc, initfracgood,initfracragged,prefix,tgtCleave)
    % initialize pool with given volume (in ul) and concentration (in nM)
    % initfracgood is the fraction of the total pool that is "good"
    % initfracragged is a two-element vector providing the fraction of the pool that is ragged on the left and right ends
    % prefix is the prefix of the library;  'W','AW','BW','-'(cleaved);  just for keeping track
    % tgtCleave is the cleavage without and with target for the "good" sequences      
      obj.volume=initvol;
      obj.nragged=obj.moles(initvol,initconc)*initfracragged;
      obj.ngood=obj.moles(initvol,initconc)*initfracgood*(1-sum(initfracragged));
      obj.nbad=obj.moles(initvol,initconc)*(1-initfracgood)*(1-sum(initfracragged));
      obj.initfracgood=initfracgood*(1-sum(initfracragged));
      obj.initngood=obj.ngood;
      obj.goodseqs=1:10000;
      obj.prefix=prefix;
      obj.tgtCleave=sort(tgtCleave);
      obj.history=[];
      fprintf('Initial prefix=%s, Target cleavage=[%.2f, %.2f], Initial good=%.3g ragged=[%.2g %.2g] bad=%.2g\n',obj.prefix,obj.tgtCleave,obj.ngood, obj.nragged, obj.nbad);
      obj.printdiv('Initial');
      obj.cumcost=0;
    end

    function k=kgood(obj)
    % Compute the mean number of copies of each good sequence
      k=length(obj.goodseqs)/length(unique(obj.goodseqs));
    end
    
    function g=divtarget(obj)
    % Compute the fraction of the diversity target we are at
      g=length(unique(obj.goodseqs))/10000 * obj.initngood;
    end
    
    function t=total(obj)
      t=obj.ngood+obj.nbad+sum(obj.nragged);
    end

    function f=fracgood(obj)
      f=obj.ngood/obj.total();
    end
    
    function c=conc(obj)
    % Concentration in nM
      c=obj.total()/(obj.volume*1e-6)/6.022e23*1e9;
    end
        
    function printdiv(obj,note)
      fprintf('%-50.50s   %2s %4.0ful %5.0fnM Total=%7.2g Ragged=[%2.0f%%,%2.0f%%] Enrich=%7.2g kgood=%8.2f goodseqs=%5.3g cumcost=$%3.0f\n',note,obj.prefix, obj.volume,obj.conc(),obj.total(),obj.nragged/obj.total()*100,obj.fracgood()/obj.initfracgood, obj.kgood(), obj.divtarget(), obj.cumcost);
      obj.history=[obj.history,struct('ngood',obj.ngood,'bad',obj.nbad,'ragged',obj.nragged,'divtarget',obj.divtarget,'note',note)];
    end

    function plothistory(obj)
      subplot(211);
      h=semilogy([obj.history.bad],'-o');
      hold on;
      h(2)=semilogy(arrayfun(@(z) z.ragged(1),obj.history),'-o');
      h(3)=semilogy(arrayfun(@(z) z.ragged(2),obj.history),'-o');
      set(gca,'XTick',1:length(obj.history));
      set(gca,'XTickLabel',{});
      legend(h,{'bad','left ragged','right ragged'});
      
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
    
    function resample(obj,note,gain)
    % Resample (with replacement) pool with given gain
      obj.ngood=obj.ngood*gain;
      obj.nbad=obj.nbad*gain;
      obj.nragged=obj.nragged*gain;
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
    
    function randchoose(obj,note,goodgain,badgain,raggedgain)
    % Choose from the pool without replacement
    % A different gain can be applied to each molecule class, each of which should be between 0 and 1
      if nargin<4
        badgain=goodgain;
        raggedgain=goodgain;
      end
      if badgain>1.0 || badgain<0.0 || goodgain>1.0 || goodgain<0.0 || any(raggedgain<0) || any(raggedgain>1)
        error('Bad gains: %f, %f, %f, %f\n', goodgain, badgain, raggedgain);
      end
      % Assume that there are exactly k copies of each sequence (TODO: does this reasonably approximate?)
      obj.ngood=obj.ngood*goodgain;
      obj.nbad=obj.nbad*badgain;
      obj.nragged=obj.nragged*raggedgain;
      obj.goodseqs=randsample(obj.goodseqs,round(goodgain*length(obj.goodseqs)),false);
      if goodgain==badgain
        obj.printdiv(sprintf('%s(gain=%.3f)',note,goodgain));
      else
        obj.printdiv(sprintf('%s(gain=%.3fG,%.3fB,%.3fR)',note,goodgain,badgain,raggedgain));
      end
    end

% Select using given cleavage rate for BAD sequences
    function Select(obj,keepCleave,cleaveRate)
    % Choose molecules from the pool that either cleave or not (depending on keepCleave)
    % The cleavage rate of the bad molecules is 'cleaveRate'.  
    % Ragged molecules never cleave and good molecules cleave with tgtCleave
      badgain=cleaveRate;
      if keepCleave
        % Compute badgain such that observed cleavage rate is obtained
        obj.prefix='-';
        obj.randchoose(sprintf('Select(Clv,clv=%.2f)',cleaveRate),obj.tgtCleave(2),badgain,0);
      else
        obj.prefix='W';
        obj.randchoose(sprintf('Select(Unclv,clv=%.2f)',cleaveRate),1-obj.tgtCleave(1),badgain,1.0);
      end
    end
    
    function T7(obj,rnaconc)
    % Transcribe the pool ending with the given RNA concentration
    % Assume molecules with ragged left end won't transcribe, but ragged right end ones are inherited
    % After this method, the pool will be referring to the RNA produced only, not including the template
      obj.cumcost=obj.cumcost+0.079*obj.volume;   % Price/ul of rx:  T7: .050 NTP .016, SuperaseIn .013
      gain=rnaconc*obj.total()/(obj.conc()*(obj.ngood+obj.nbad+obj.nragged(2)));
      obj.nragged(1)=0;
      obj.resample('T7',gain);
    end

    function RT(obj, efficiency,primerRaggedFrac)
      obj.cumcost=obj.cumcost+.236*obj.volume;
      newRagged=(obj.ngood+obj.nbad)*efficiency*primerRaggedFrac;
      gain=efficiency*(1-primerRaggedFrac);
      obj.nragged(2)=newRagged/gain;
      obj.randchoose('RT',gain);
    end
    
    function ligate(obj, prefix, efficiency)
      obj.prefix=prefix;
      obj.randchoose(sprintf('ligate(%s)',prefix),efficiency);
      obj.cumcost=obj.cumcost+.0002*obj.volume;
    end
    
    function gain=PCR(obj,ncycles,primersRaggedFrac)
    % PCR amplify the pool to the given final volume and concentration
    % Assumes perfect amplication and uniform copying of all non-ragged input molecules
    % The resulting pool will inherit the ragged fraction of the primers + the original unamplified ragged ones
      S1=struct('nbad',obj.nbad,'ngood',obj.ngood,'nragged',[obj.nragged(1),0]);
      S2=struct('nbad',obj.nbad,'ngood',obj.ngood,'nragged',[obj.nragged(2),0]);
      for i=1:ncycles
        % Update
        fprintf('S1=%.2g %.2g [%.2g %.2g] S2=%.2g %.2g [%.2g %.2g]\n', S1.nbad, S1.ngood, S1.nragged, S2.nbad, S2.ngood, S2.nragged);
        S1new=struct('nbad',S1.nbad+S2.nbad*(1-primersRaggedFrac(1)),...
                     'ngood',S1.ngood+S2.ngood*(1-primersRaggedFrac(1)),...
                     'nragged',S1.nragged+[primersRaggedFrac(1)*(S2.ngood+S2.nbad),S2.nragged(1)]);
        S2new=struct('nbad',S2.nbad+S1.nbad*(1-primersRaggedFrac(1)),...
                     'ngood',S2.ngood+S1.ngood*(1-primersRaggedFrac(1)),...
                     'nragged',S2.nragged+[primersRaggedFrac(1)*(S1.ngood+S1.nbad),S1.nragged(1)]);
        S1=S1new;
        S2=S2new;
        obj.goodseqs=[obj.goodseqs,randsample(obj.goodseqs,round(length(obj.goodseqs)*(1-primersRaggedFrac(1))),false)];
      end
      fprintf('S1=%.2g %.2g [%.2g %.2g] S2=%.2g %.2g [%.2g %.2g]\n', S1.nbad, S1.ngood, S1.nragged, S2.nbad, S2.ngood, S2.nragged);
      % No change to number of good sequences since sampling is uniform
      oldtotal=obj.total;
      obj.ngood=min(S1.ngood,S2.ngood);
      obj.nbad=min(S1.nbad,S2.nbad);
      obj.nragged=[max(S1.nragged(1),S2.nragged(2)),max(S1.nragged(2),S2.nragged(1))];
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
  end
end
