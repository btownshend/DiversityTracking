% DivTrack - track diversity of a series of rounds of TRP

classdef DivTrack < matlab.mixin.Copyable
  properties
    volume;
    conc;
    kgood;
    kbad;
    fracgood;
    initfracgood;
    ngoodinitial;
    tgtCleave;
    prefix;
    initdiversity;
  end

  methods(Static)
    function m=moles(vol,conc)
      m=vol*1e-6*conc*1e-9*6.022e23;
    end
  end
  
  methods
    function obj=DivTrack(initvol, initconc, initfracgood,prefix,tgtCleave)
      obj.volume=initvol;
      obj.conc=initconc;
      obj.kgood=1;
      obj.kbad=1;
      obj.initfracgood=initfracgood;
      obj.fracgood=initfracgood;
      obj.ngoodinitial=obj.moles(obj.volume,obj.conc)*obj.fracgood;
      obj.prefix=prefix;
      obj.tgtCleave=sort(tgtCleave);
      obj.initdiversity=obj.moles(initvol,initconc);
      fprintf('Initial prefix=%s, Target cleavage=[%.2f, %.2f], Initial diversity=%.2g Fraction good=%.2g\n',obj.prefix,obj.tgtCleave,obj.initdiversity,obj.fracgood);
      obj.printdiv('Initial');
    end

    function printdiv(obj,note)
      fprintf('%-40.40s   %c %4.0ful %4.0fnM Div=%7.2g Enrich=%7.2g (%7.2g good kept) kgood=%8.2f kbad=%8.2f fracgood=%.1g\n',note,obj.prefix, obj.volume,obj.conc,obj.moles(obj.volume,obj.conc)/obj.kgood,obj.fracgood/obj.initfracgood, obj.moles(obj.volume,obj.conc)/obj.kgood*obj.fracgood/obj.ngoodinitial,obj.kgood, obj.kbad, obj.fracgood);
    end

    function resample(obj,note,goodgain,badgain)
      if nargin<4
        badgain=goodgain;
      end
      obj.kgood=obj.kgood*goodgain/(1-exp(-obj.kgood*goodgain));
      obj.kbad=obj.kbad*badgain/(1-exp(-obj.kbad*badgain));
      obj.fracgood=obj.fracgood*goodgain/(obj.fracgood*goodgain+(1-obj.fracgood)*badgain);
      obj.conc=obj.conc*(goodgain*obj.fracgood + badgain*(1-obj.fracgood));
      if goodgain==badgain
        obj.printdiv(sprintf('%s(gain=%.3f)',note,goodgain));
      else
        obj.printdiv(sprintf('%s(gain=%.3fG,%.3fB)',note,goodgain,badgain));
      end
    end

    function randchoose(obj,note,goodgain,badgain)
      if nargin<4
        badgain=goodgain;
      end
      if badgain>1.0 || badgain<0.0 || goodgain>1.0 || goodgain<0.0
        error('Bad gains: %f, %f\n', goodgain, badgain);
      end
      % Assume that there are exactly k copies of each sequence (TODO: does this reasonably approximate?)
      pmissedgood=(1-goodgain)^obj.kgood;
      obj.kgood=obj.kgood*goodgain/(1-pmissedgood);
      pmissedbad=(1-badgain)^obj.kbad;
      obj.kbad=obj.kbad*badgain/(1-pmissedbad);
      obj.fracgood=obj.fracgood*goodgain/(obj.fracgood*goodgain+(1-obj.fracgood)*badgain);
      obj.conc=obj.conc*(goodgain*obj.fracgood + badgain*(1-obj.fracgood));
      if goodgain==badgain
        obj.printdiv(sprintf('%s(gain=%.3f)',note,goodgain));
      else
        obj.printdiv(sprintf('%s(gain=%.3fG,%.3fB)',note,goodgain,badgain));
      end
    end

    function Select(obj,keepCleave,cleaveRate)
      if keepCleave
        obj.randchoose(sprintf('Select(Clv,clv=%.2f)',cleaveRate),obj.tgtCleave(2),cleaveRate);
        if obj.prefix=='A'
          obj.prefix='B';
        else
          obj.prefix='A';
        end
      else
        obj.randchoose(sprintf('Select(Unclv,clv=%.2f)',cleaveRate),1-obj.tgtCleave(1),1-cleaveRate);
      end
    end
    
    function T7(obj,rnaconc)
      obj.resample('T7',rnaconc/obj.conc);
    end

    function gain=PCR(obj,finalvol,finalconc)
      obj.conc=obj.conc/(finalvol/obj.volume);
      obj.volume=finalvol;
      gain=obj.moles(finalvol,finalconc)/obj.moles(obj.volume,obj.conc);
      obj.resample('PCR',gain);
    end

    function measured(obj,note,volume,conc)
      obj.conc=obj.conc/(volume/obj.volume);
      obj.volume=volume;
      gain=obj.moles(volume,conc)/obj.moles(obj.volume,obj.conc);
      obj.randchoose(note,gain);
    end
  end
end
