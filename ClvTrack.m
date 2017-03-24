% ClvTrack - track precleavage selection

classdef ClvTrack < matlab.mixin.Copyable
  properties
    volume; % Volume in ul
    ndna;   % ndna(type) - see constants below for indexing
    cumcost;
    trackfrac;	% Number of tracked molecules per total molecules
    cdist;	% array of ClvDistributions indexed by TYPE
    rnaLength;  % Uncleaved RNA length in nt
    mw;		% Molecular weight of each type
    active;	% Currently active stage (i.e. which is the relevant cdist)
    id;		% ID of current sample
  end

  properties (Constant)
    TYPE_W=1;
    TYPE_NP=2;
    TYPE_T7Wf=3;
    TYPE_T7Wr=4;
    TYPE_CIRC=5;
    TYPE_URNA=6;
    TYPE_CRNA=7;
    TYPE_NTYPES=7;
  end
  
  methods(Static)
    function m=moles(vol,conc)
      m=vol*1e-6*conc*1e-9*6.022e23;
    end

    function dsdna=MWdsDNA(len)
      dsdna=len*607.4+157.9;
      ssdna=len*303.7+79;
      rdna=len*320.5+159.0;
    end
    function ssdna=MWssDNA(len)
      ssdna=len*303.7+79;
    end
    function rna=MWRNA(len)
      rna=len*320.5+159.0;
    end
    
    function ngul=conc2ngul(conc,mw)
      ngul=conc*1e-9*mw*1000;
    end
   
  end
  
  
  methods
    function obj=ClvTrack(initvol, initconc, cleavages, id, note, trackfrac)
    % initialize pool with given volume (in ul) and concentration (in nM)
    % cleavages are a vector sample of observed cleavages (actually fitness over a single round)
      obj.volume=initvol;
      obj.ndna=zeros(obj.TYPE_NTYPES,1);
      obj.ndna(obj.TYPE_T7Wr)=obj.moles(initvol,initconc);
      if nargin>=6
        obj.trackfrac=trackfrac;
      else
        obj.trackfrac=1000/sum(obj.ndna);
      end
      obj.cdist=cell(1,obj.TYPE_NTYPES);
      cleavelist=randsample(cleavages,ceil(obj.trackfrac*max(obj.ndna)),true);
      for i=1:obj.TYPE_NTYPES
        obj.cdist{i}=ClvDistribution(cleavelist,round(obj.trackfrac*obj.ndna(i)));
      end
      obj.active=obj.TYPE_T7Wr;
      if nargin>=4
        obj.id=id;
      else
        obj.id='';
      end
      if nargin<5
        note='Initial';
      end
      obj.printdiv(note,true);
      obj.cumcost=initvol*initconc/1e6*235;   % Template cost
      obj.rnaLength=108;
      extLength=71;   % Length added to cDNA beyond RNA end due to RT primer
      obj.mw=zeros(1,obj.TYPE_NTYPES);
      obj.mw(obj.TYPE_W)=obj.MWssDNA(obj.rnaLength+extLength);
      obj.mw(obj.TYPE_NP)=obj.MWssDNA(obj.rnaLength+extLength-20);
      obj.mw(obj.TYPE_T7Wr)=obj.MWdsDNA(obj.rnaLength+21)/2;
      obj.mw(obj.TYPE_T7Wf)=obj.MWdsDNA(obj.rnaLength+21)/2;
      obj.mw(obj.TYPE_CIRC)=obj.mw(obj.TYPE_NP);
      obj.mw(obj.TYPE_URNA)=obj.MWRNA(obj.rnaLength);
      obj.mw(obj.TYPE_CRNA)=obj.mw(obj.TYPE_URNA);   % Counts the 5' cleavage product too
    end

    function c=deepcopy(obj)
      c=obj.copy();
      for i=1:length(c.cdist)
        c.cdist{i}=c.cdist{i}.copy();
      end
    end
    
    function changetrackfrac(obj,newfrac)
    % Reduce the tracking fraction to match o2
      assert(false);   % TODO: UPDATE
      fprintf('Updating tracking fraction of %s from %g to %g (by %g)\n', obj.id, obj.trackfrac, newfrac, newfrac/obj.trackfrac);
      minid=inf;
      maxid=0;
      for i=1:length(obj.goodseqs)
        minid=min([minid,obj.goodseqs{i}]);
        maxid=max([maxid,obj.goodseqs{i}]);
      end
      while newfrac>obj.trackfrac
        % replicate each entry
        offset=maxid-minid+1;
        for i=1:length(obj.goodseqs)
          obj.goodseqs{i}=[obj.goodseqs{i},obj.goodseqs{i}+offset];
        end
        maxid=maxid+offset;
        obj.trackfrac=obj.trackfrac*2;
      end
      if (newfrac<obj.trackfrac)
        % Reduce maxid
        maxid=round((maxid-minid+1)*newfrac/obj.trackfrac+minid-1);
        for i=1:length(obj.goodseqs)
          obj.goodseqs{i}=obj.goodseqs{i}(obj.goodseqs{i}<=maxid);
        end
      end
      obj.trackfrac=newfrac;
      obj.printdiv('After updating trackfrac');
    end
  
    function mix(obj,o2)
    % Add another tracked set to this one
      if o2.trackfrac<obj.trackfrac
        o2.changetrackfrac(obj.trackfrac);
      elseif obj.trackfrac<o2.trackfrac
        obj.changetrackfrac(o2.trackfrac);
      end

      obj.cumcost=obj.cumcost+o2.cumcost;
      
      offset=ceil(max(cellfun(@(z) max([z.seqs,0]), obj.cdist))/1000)*1000;
      for i=1:length(obj.cdist)
        %fprintf('%d/%d + %d/%d ->', length(obj.cdist{i}.seqs),length(obj.cdist{i}.cleavage), length(o2.cdist{i}.seqs),length(o2.cdist{i}.cleavage));
        obj.cdist{i}.seqs=[obj.cdist{i}.seqs,o2.cdist{i}.seqs+offset];
        obj.cdist{i}.cleavage((1:length(o2.cdist{i}.cleavage))+offset)=o2.cdist{i}.cleavage;
        %fprintf('%d/%d\n', length(obj.cdist{i}.seqs),length(obj.cdist{i}.cleavage));
        obj.ndna(i)=obj.cdist{i}.nseq/obj.trackfrac;
      end
      obj.volume=obj.volume+o2.volume;
      obj.printdiv(sprintf('Mix(%s)',o2.id));
    end

    function k=avgcopies(obj)
    % Compute the mean number of copies of each good sequence
      [nu,n]=obj.cdist{obj.active}.nunique();
      k=n/nu;
    end
    
    function t=total(obj)
      t=sum(obj.ndna(:));
    end

    function c=conc(obj,typ)
    % Concentration in nM
      if nargin<2 || isempty(typ)
        typ=1:size(obj.ndna,1);
      end
      c=sum(obj.ndna(typ))/(obj.volume*1e-6)/6.022e23*1e9;
    end
        
    function printdiv(obj,note, heading)
      if nargin<3
        heading=false;
      end
      concs='';
      marker=' *';
      for i=1:size(obj.ndna,1)
        concs=[concs,sprintf('%6.1f%c',obj.conc(i),marker((i==obj.active)+1))];
      end
      highcleave=.85;
      if heading
        fprintf('\n%-7.7s%-50.50s    %5s %6s %6s %6s %6s %6s %6s %6s      Total  Clv  kCopy  pmoled  k>%2.0f%%   p>%2.0f%%  NTrk Cost\n','ID','Desc','Volume','W','NoPre','T7Wf','T7Wr','Circ','URNA','CRNA',100*highcleave,100*highcleave);
        fprintf('%s\n',repmat('-',173,1));
      end
      pmolediv=obj.ndna(obj.active)/obj.avgcopies()*1e12/6.022e23;
      frachigh=obj.cdist{obj.active}.frachigh(highcleave);
      khigh=obj.cdist{obj.active}.khigh(highcleave);
      fprintf('%-7.7s%-50.50s   %5.0ful %snM %7.2g %4.1f %6.2f %7.2f %6.2f %7.2f %5.0f $%3.0f\n',...
              obj.id,note,obj.volume,concs,obj.total(), 100*obj.cdist{obj.active}.meanCleavage(), ...
              obj.avgcopies(), pmolediv, khigh, pmolediv*frachigh, obj.cdist{obj.active}.nunique(), obj.cumcost);
      for i=1:length(obj.cdist)
        if abs(obj.ndna(i)*obj.trackfrac-obj.cdist{i}.nseq())>1
          fprintf('ndna(%d,1)*trackfrac=%.0f, length(cdist)=%d\n', i, obj.ndna(i)*obj.trackfrac, obj.cdist{i}.nseq());
          keyboard;
        end
      end
    end

    function resample(obj,note,gain)
    % Resample (with replacement) pool with given gain
      obj.ndna=obj.ndna*gain;
      for i=1:obj.TYPE_NTYPES
        obj.cdist{i}.resample(round(obj.ndna(i)*obj.trackfrac),true);
      end
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
    
    function randchoose(obj,note,gain)
    % Choose from the pool without replacement
      if gain>1.0 || gain<0.0
        error('Bad gain: %f\n', gain);
      end
      % Assume that there are exactly k copies of each sequence (TODO: does this reasonably approximate?)
      obj.ndna=obj.ndna*gain;
      for i=1:obj.TYPE_NTYPES
        nsamp=round(obj.ndna(i)*obj.trackfrac);
        obj.cdist{i}.resample(nsamp,false);
      end
      obj.printdiv(sprintf('%s(gain=%.3f)',note,gain));
    end

    function T7(obj,rnaconc)
    % Transcribe the pool ending with the given RNA concentration
    % After this method, the pool will be referring to the RNA produced only, not including the template
      obj.cumcost=obj.cumcost+0.079*obj.volume;   % Price/ul of rx:  T7: .050 NTP .016, SuperaseIn .013
      gain=rnaconc/obj.conc(obj.TYPE_T7Wr);
      rna=obj.cdist{obj.TYPE_T7Wr}.copy();
      rna.resample(round(length(rna.seqs)*gain),true);
      [obj.cdist{obj.TYPE_URNA},obj.cdist{obj.TYPE_CRNA}]=rna.splitBasedOnCleavage();
      obj.ndna(obj.TYPE_URNA)=obj.cdist{obj.TYPE_URNA}.nseq()/obj.trackfrac;
      obj.ndna(obj.TYPE_CRNA)=obj.cdist{obj.TYPE_CRNA}.nseq()/obj.trackfrac;
      obj.active=obj.TYPE_CRNA;
      obj.printdiv(sprintf('T7(gain=%.3f,cleavage=%.0f%%)',gain, 100*obj.ndna(obj.TYPE_CRNA)/sum(obj.ndna([obj.TYPE_URNA,obj.TYPE_CRNA]))));
    end

    function RT(obj, efficiency)
    % Reverse transcription
      stopconc=2000;
      obj.cumcost=obj.cumcost+.236*obj.volume;
      rnaconc=obj.conc(obj.TYPE_URNA)+obj.conc(obj.TYPE_CRNA);
      if rnaconc*efficiency>stopconc
        %efficiency=stopconc/rnaconc;
        fprintf('Should be limiting RT efficiency to %.2f due to limited stop oligo\n', stopconc/rnaconc);
      end
      % Apply the gain
      nsamp=round(obj.ndna(obj.TYPE_URNA)*efficiency*obj.trackfrac);
      obj.cdist{obj.TYPE_W}.addsample(obj.cdist{obj.TYPE_URNA},nsamp,false);
      obj.cdist{obj.TYPE_URNA}.clear();
      obj.ndna(obj.TYPE_W)=obj.ndna(obj.TYPE_W)+obj.ndna(obj.TYPE_URNA)*efficiency;
      obj.ndna(obj.TYPE_URNA)=0;  % Degrade RNA
      
      nsamp=round(obj.ndna(obj.TYPE_CRNA)*efficiency*obj.trackfrac);
      obj.cdist{obj.TYPE_NP}.addsample(obj.cdist{obj.TYPE_CRNA},nsamp,false);
      obj.ndna(obj.TYPE_NP)=obj.ndna(obj.TYPE_NP)+obj.ndna(obj.TYPE_CRNA)*efficiency;
      obj.cdist{obj.TYPE_CRNA}.clear();
      obj.ndna(obj.TYPE_CRNA)=0;  % Degrade RNA

      obj.active=obj.TYPE_NP;
      obj.printdiv(sprintf('RT(eff=%.3f)',efficiency));
    end
    
    function ligate(obj, efficiency)
      nligated=obj.ndna(obj.TYPE_NP)*efficiency;
      [c1,obj.cdist{obj.TYPE_NP}]=obj.cdist{obj.TYPE_NP}.splitRandom(efficiency);
      obj.cdist{obj.TYPE_CIRC}.addsample(c1,length(c1.seqs),false);
      obj.ndna(obj.TYPE_NP)=obj.cdist{obj.TYPE_NP}.nseq()/obj.trackfrac;
      obj.ndna(obj.TYPE_CIRC)=obj.cdist{obj.TYPE_CIRC}.nseq()/obj.trackfrac;
      obj.active=obj.TYPE_CIRC;
      obj.printdiv(sprintf('Ligate(eff=%.3f)',efficiency));
      obj.cumcost=obj.cumcost+.0002*obj.volume;
    end
    
    function exo(obj,residual)
    % Exo digestion -- remove all DNA but circ
      if nargin<2
        residual=0;
      end
      remove=[obj.TYPE_W,obj.TYPE_NP,obj.TYPE_T7Wr,obj.TYPE_T7Wf];
      % Some residual of T7W will appear as W prefix
      obj.ndna(obj.TYPE_W)=obj.ndna(obj.TYPE_T7Wr)*residual;
      obj.ndna([obj.TYPE_T7Wr,obj.TYPE_T7Wf,obj.TYPE_NP])=0;
      obj.cdist{obj.TYPE_W}.clear();
      obj.cdist{obj.TYPE_W}.addsample(obj.cdist{obj.TYPE_T7Wr},round(obj.ndna(obj.TYPE_W)*obj.trackfrac),false);
      obj.cdist{obj.TYPE_T7Wr}.clear();
      obj.cdist{obj.TYPE_T7Wf}.clear();
      obj.cdist{obj.TYPE_NP}.clear();
      obj.printdiv(sprintf('Exo(resid=%.2f)',residual));
    end
    
    function printmeasure(obj,note,obs,err)
      allconc=nan(size(obs));
      for i=1:length(obs)
        if isfinite(obs(i))
          allconc(i)=obj.conc(i);
        end
      end
      allrel=obs./allconc;
      allrelerr=allrel;
      allrelerr(allrelerr<1)=1./allrelerr(allrelerr<1);
      rel=nanmean(allrel);
      fprintf('%-7.7s%50.50s           %snM ','',note,sprintf('%6.1f ',obs));
      for i=1:length(err)
        rel=err(i);
        if isfinite(rel)
          str=repmat('*',1,round(rel*10));
          if length(str)>10
            str(10)='|';
          end
          if length(str)>20
            str(20)='|';
          end
          str=['|',str];
        else
          str='';
        end
        fprintf('%4.2f %-25.25s ',rel, str);
      end
      fprintf('\n');
    end
    
    function qpcr(obj,note,t7,w,m)
    % Note qPCR measurement (in nM)
      if nargin<4
        w=nan;
      end
      if nargin<5
        m=nan;
      end
      obs=nan(size(obj.ndna,1),1);
      rel=obj.ndna([obj.TYPE_CIRC,obj.TYPE_T7Wr,obj.TYPE_T7Wf])/sum(obj.ndna([obj.TYPE_CIRC,obj.TYPE_T7Wr,obj.TYPE_T7Wf]));
      obs(obj.TYPE_CIRC)=t7*rel(1);
      obs(obj.TYPE_T7Wr)=t7*rel(2);
      obs(obj.TYPE_T7Wf)=t7*rel(3);
      if isfinite(m)
        rel=obj.ndna([obj.TYPE_W,obj.TYPE_NP]);
        rel=rel/sum(rel);
        obs(obj.TYPE_W)=(m-t7)*rel(1);
        obs(obj.TYPE_NP)=(m-t7)*rel(2);
      end
      
      err=t7/nansum(obj.ndna([obj.TYPE_CIRC,obj.TYPE_T7Wr,obj.TYPE_T7Wf]));
      if isfinite(w)
        %err(2)=w/nansum(nansum(obj.ndna([obj.TYPE_CIRC,obj.TYPE_T7Wr,obj.TYPE_T7Wf,obj.TYPE_W])));
        % W doesn't reflect true W since stop oligo blocks extension
        err(2)=0;
      end
      if isfinite(m)
        err(3)=m/nansum(obj.ndna([obj.TYPE_CIRC,obj.TYPE_T7Wr,obj.TYPE_T7Wf,obj.TYPE_W,obj.TYPE_NP]));
      end
      err=err*1e-9*obj.volume*1e-6*6.022e23;
      obj.printmeasure(sprintf('qPCR(%.1f,%.1f,%.1f) %s',t7,w,m,note),obs,err);
    end
    
    function qubitdna(obj,note,ngul)
    % Qubit HS DNA concentration as given
      ssFactor=1.0;     % Assume ssDNA reads at ssFactor of dsDNA
      factors=[ssFactor ssFactor 1 1 ssFactor 0 0];
      allconc=[];
      for i=1:size(obj.ndna,1)
        allconc(i)=obj.conc(i);
      end
      expect=allconc.*factors.*obj.mw*1e-9*1000;
      scale=ngul/nansum(expect);
      obs=expect*scale./factors./obj.mw/1e-9/1000;
      obs(factors==0)=nan;
      obj.printmeasure(sprintf('qubit HS DNA(%.1f ng/ul) %s',ngul,note),obs,scale);
    end

    function qubitrna(obj,note,ngul)
    % Qubit HS RNA concentration as given
      factors=[0 0 0 0 0 1 1];
      allconc=[];
      for i=1:size(obj.ndna,1)
        allconc(i)=obj.conc(i);
      end
      expect=allconc.*factors.*obj.mw*1e-9*1000;
      scale=ngul/nansum(expect);
      obs=expect*scale./factors./obj.mw/1e-9/1000;
      obs(factors==0)=nan;
      obj.printmeasure(sprintf('qubit RNA(%.1f ng/ul) %s',ngul,note),obs,scale);
    end

    function nanodrop(obj,note,ngul,r1,r2)
      ratios='';
      if nargin>=4
        ratios=[ratios,sprintf(' %.2f',r1)];
      end
      if nargin>=5
        ratios=[ratios,sprintf(' %.2f',r2)];
      end
      a260=ngul/50;    % Assume was set for DNA-50
      factors=[1/33 1/33 1/50 1/50 1/33 1/40 1/40];
      allconc=[];
      for i=1:size(obj.ndna,1)
        allconc(i)=obj.conc(i);
      end
      expect=allconc.*factors.*obj.mw*1e-9*1000;
      scale=a260/nansum(expect);
      obs=expect*scale./factors./obj.mw/1e-9/1000;
      obs(factors==0)=nan;
      obj.printmeasure(sprintf('Nanodrop(%.1f ng/ul%s) %s',ngul,ratios,note),obs,scale);
    end
      
    function gain=PCR(obj,ncycles,maxconc)
    % PCR amplify the pool to the given final volume and concentration
    % Assumes perfect amplication and uniform copying of all input molecules
      if nargin<3
        maxconc=175;
      end
      maxn=maxconc*1e-9*6.022e23*obj.volume*1e-6;
      oldtotal=obj.total();
      for i=1:ncycles
        n1=obj.ndna(obj.TYPE_T7Wf);
        if n1+obj.ndna(obj.TYPE_T7Wr)>maxn
          n1=maxn-obj.ndna(obj.TYPE_T7Wr);
        end
        if n1>0
          nr=obj.ndna(obj.TYPE_T7Wr)+n1;
          tmpdist=obj.cdist{obj.TYPE_T7Wr}.copy();
          tmpdist.addsample(obj.cdist{obj.TYPE_T7Wf},round(n1*obj.trackfrac),false);
        end
        n2f=obj.ndna(obj.TYPE_T7Wr);
        n2c=obj.ndna(obj.TYPE_CIRC);
        if n2f+n2c+obj.ndna(obj.TYPE_T7Wf)>maxn
          g2=(maxn-obj.ndna(obj.TYPE_T7Wf))/(n2f+n2c);
          n2f=n2f*g2;
          n2c=n2c*g2;
        end
        if n2f+n2c>0
          obj.cdist{obj.TYPE_T7Wf}.addsample(obj.cdist{obj.TYPE_T7Wr},round(n2f*obj.trackfrac),false);
          obj.cdist{obj.TYPE_T7Wf}.addsample(obj.cdist{obj.TYPE_CIRC},round(n2c*obj.trackfrac),false);
          obj.ndna(obj.TYPE_T7Wf)=obj.cdist{obj.TYPE_T7Wf}.nseq()/obj.trackfrac;
        end
        if n1>0
          obj.cdist{obj.TYPE_T7Wr}=tmpdist;
          obj.ndna(obj.TYPE_T7Wr)=tmpdist.nseq()/obj.trackfrac;
        end
      end
      gain=obj.total()/oldtotal;
      obj.cumcost=obj.cumcost+obj.volume*263/250/50;  % Kapa is $263/250U, uses 1U/50ul
      obj.active=obj.TYPE_T7Wr;
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

    function cleanup(obj,note,volume,efficiency)
      obj.volume=volume;
      obj.randchoose(note,efficiency);
    end
  end
end
