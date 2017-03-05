% TRPTrack - track diversity of a series of rounds of TRP

classdef TRPTrack < matlab.mixin.Copyable
  properties
    volume; % Volume in ul
    ndna;   % ndna(type,goodbad) - type 1=W, 2=no prefix, 3=T7W, 4=circ, 5-W_RNA, 6-clvd RNA; goodbad 1=good, 2=bad
    history;
    cumcost;
    initfracgood;  % initfracgood is the initial fraction of the total pool that is "good"
    trackfrac;	   % Number of tracked molecules per good molecule (i.e. length(goodseqs)/ndna(GOOD) )
    goodseqs;  	   % goodseqs{type} - Each value in this vector represents a unique sequence, keeps empirical distribution
    rnaLength;  % Uncleaved RNA length in nt
    mw;		% Molecular weight of each type
    active;	% Currently active stage (i.e. which is the relevant goodseqs)
    id;		% ID of current sample
  end

  properties (Constant)
    TYPE_W=1;
    TYPE_NP=2;
    TYPE_T7W=3;
    TYPE_CIRC=4;
    TYPE_URNA=5;
    TYPE_CRNA=6;
    GOOD=1;
    BAD=2;
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
    function obj=TRPTrack(initvol, initconc, fracgood, id, note, trackfrac)
    % initialize pool with given volume (in ul) and concentration (in nM)
    % poolFitness is an empirical distribution of the cleavage of the pool, each value is equally probable
      obj.volume=initvol;
      obj.ndna=zeros(6,2);
      nmolecules=obj.moles(initvol,initconc);
      obj.ndna(obj.TYPE_T7W,obj.GOOD)=nmolecules*fracgood;
      obj.ndna(obj.TYPE_T7W,obj.BAD)=nmolecules*(1-fracgood);
      obj.initfracgood=obj.fracgood();    
      if nargin>=6
        obj.trackfrac=trackfrac;
      else
        obj.trackfrac=1000/sum(obj.ndna(:,obj.GOOD));
      end
      obj.goodseqs=cell(1,6);
      icnt=1;
      for i=1:size(obj.ndna,1)
        n=round(obj.ndna(i,obj.GOOD)*obj.trackfrac);
        obj.goodseqs{i}=icnt:icnt+n-1;
        icnt=icnt+n;
      end
      obj.active=obj.TYPE_T7W;
      obj.history=[];
      if nargin>=4
        obj.id=id;
      else
        obj.id='';
      end
      %fprintf('Initial good=%.3g bad=%.2g\n',obj.ngood(), obj.nbad());
      if nargin<5
        note='Initial';
      end
      obj.printdiv(note,true);
      obj.cumcost=initvol*initconc/1e6*235;   % Template cost
      obj.rnaLength=108;
      extLength=71;   % Length added to cDNA beyond RNA end due to RT primer
      obj.mw=zeros(1,6);
      obj.mw(obj.TYPE_W)=obj.MWssDNA(obj.rnaLength+extLength);
      obj.mw(obj.TYPE_NP)=obj.MWssDNA(obj.rnaLength+extLength-20);
      obj.mw(obj.TYPE_T7W)=obj.MWdsDNA(obj.rnaLength+21);
      obj.mw(obj.TYPE_CIRC)=obj.mw(obj.TYPE_NP);
      obj.mw(obj.TYPE_URNA)=obj.MWRNA(obj.rnaLength);
      obj.mw(obj.TYPE_CRNA)=obj.mw(obj.TYPE_URNA);   % Counts the 5' cleavage product too
    end

    function changetrackfrac(obj,newfrac)
    % Reduce the tracking fraction to match o2
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

      obj.ndna=obj.ndna+o2.ndna;
      obj.cumcost=obj.cumcost+o2.cumcost;
      obj.initfracgood=(obj.volume*obj.initfracgood+o2.volume*o2.initfracgood)/(obj.volume+o2.volume);
      
      offset=ceil(max(cellfun(@(z) max([z,0]), obj.goodseqs))/1000)*1000;
      for i=1:length(obj.goodseqs)
        obj.goodseqs{i}=[obj.goodseqs{i},o2.goodseqs{i}+offset];
      end
      obj.volume=obj.volume+o2.volume;
      obj.printdiv(sprintf('Mix(%s)',o2.id));
    end
    
    function [nu,n]=nunique(obj)
    % Compute the mean number of copies of each good sequence
      nu=length(unique(obj.goodseqs{obj.active}));
      n=length(obj.goodseqs{obj.active});
    end

    function k=kgood(obj)
    % Compute the mean number of copies of each good sequence
      [nu,n]=obj.nunique();
      k=n/nu;
    end
    
    function n=ngood(obj)
      n=sum(obj.ndna(:,obj.GOOD));
    end

    function n=nbad(obj)
      n=sum(obj.ndna(:,obj.BAD));
    end
    
    function t=total(obj)
      t=sum(obj.ndna(:));
    end

    function f=fracgood(obj)
      f=obj.ngood()/obj.total();
    end
    
    function c=conc(obj,typ,goodbad)
    % Concentration in nM
      if nargin<2 || isempty(typ)
        typ=1:size(obj.ndna,1);
      end
      if nargin<3 || isempty(goodbad)
        goodbad=1:2;
      end
      c=sum(sum(obj.ndna(typ,goodbad)))/(obj.volume*1e-6)/6.022e23*1e9;
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
      if heading
        fprintf('\n%-7.7s%-50.50s    %5s %6s %6s %6s %6s %6s %6s      Total  Enrich    kGood NGood Cost\n','ID','Desc','Volume','W','NoPre','T7W','Circ','URNA','CRNA');
        fprintf('%s\n',repmat('-',148,1));
      end
      fprintf('%-7.7s%-50.50s   %5.0ful %snM %7.2g %7.2g %8.2f %5.0f $%3.0f\n',obj.id,note,obj.volume,concs,obj.total(),obj.fracgood()/obj.initfracgood, obj.kgood(), obj.nunique(), obj.cumcost);
      obj.history=[obj.history,struct('ngood',obj.ngood(),'bad',obj.nbad(),'nunique',obj.nunique,'note',note)];
      for i=1:length(obj.goodseqs)
        if abs(obj.ndna(i,1)*obj.trackfrac-length(obj.goodseqs{i}))>1  && abs(obj.ndna(i,1)*obj.trackfrac/length(obj.goodseqs{i})-1)>0.01
          fprintf('ndna(%d,1)*trackfrac=%.0f, length(goodseqs)=%d\n', i, obj.ndna(i,1)*obj.trackfrac, length(obj.goodseqs{i}));
          keyboard;
        end
      end
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
      h(2)=semilogy([obj.history.nunique],'-o');
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
      dist=dist/sum(dist)*obj.nunique();
      setfig(ti);clf;
      plot(n(n>0),dist(n>0),'o-');
      xlabel('Number of copies of sequence');
      ylabel('Number of sequences');
      title(sprintf('%s: Total of %.2f sequences with mean of %.2f copies',ti, obj.nunique(), obj.kgood()));
    end
    
    function resample(obj,note,gain)
    % Resample (with replacement) pool with given gain
      obj.ndna=obj.ndna*gain;
      for i=1:length(obj.goodseqs)
        obj.goodseqs{i}=randsample(obj.goodseqs{i},round(gain*length(obj.goodseqs{i})),true);
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
      obj.ndna(:,obj.GOOD)=obj.ndna(:,obj.GOOD)*goodgain;
      obj.ndna(:,obj.BAD)=obj.ndna(:,obj.BAD)*badgain;
      for i=1:length(obj.goodseqs)
        nsamp=round(obj.ndna(i,obj.GOOD)*obj.trackfrac);
        assert(nsamp<=length(obj.goodseqs{i}));
        if nsamp==0
          obj.goodseqs{i}=[];
        elseif length(obj.goodseqs{i})>1
          obj.goodseqs{i}=randsample(obj.goodseqs{i},nsamp,false);
        end
      end
      if goodgain==badgain
        obj.printdiv(sprintf('%s(gain=%.3f)',note,goodgain));
      else
        obj.printdiv(sprintf('%s(gain=%.3fG,%.3fB)',note,goodgain,badgain));
      end
    end

    function T7(obj,rnaconc,goodcleavage,badcleavage)
    % Transcribe the pool ending with the given RNA concentration
    % After this method, the pool will be referring to the RNA produced only, not including the template
      obj.cumcost=obj.cumcost+0.079*obj.volume;   % Price/ul of rx:  T7: .050 NTP .016, SuperaseIn .013
      gain=rnaconc/obj.conc(obj.TYPE_T7W);
      obj.ndna(obj.TYPE_URNA,obj.GOOD)=obj.ndna(obj.TYPE_T7W,obj.GOOD)*gain*(1-goodcleavage);
      obj.ndna(obj.TYPE_CRNA,obj.GOOD)=obj.ndna(obj.TYPE_T7W,obj.GOOD)*gain*goodcleavage;
      bulkcleavage=(badcleavage*obj.ndna(obj.TYPE_T7W,obj.BAD)+goodcleavage*obj.ndna(obj.TYPE_T7W,obj.GOOD))/sum(obj.ndna(obj.TYPE_T7W,:));
      %fprintf('good=%g, bulk=%g, bad=%g\n',goodcleavage, bulkcleavage, badcleavage);
      obj.ndna(obj.TYPE_URNA,obj.BAD)=obj.ndna(obj.TYPE_T7W,obj.BAD)*gain*(1-badcleavage);
      obj.ndna(obj.TYPE_CRNA,obj.BAD)=obj.ndna(obj.TYPE_T7W,obj.BAD)*gain*badcleavage;
      % Resample (with replacement) pool with given gain
      obj.goodseqs{obj.TYPE_URNA}=[obj.goodseqs{obj.TYPE_URNA},randsample(obj.goodseqs{obj.TYPE_T7W},round((1-goodcleavage)*gain*length(obj.goodseqs{obj.TYPE_T7W})),true)];
      obj.goodseqs{obj.TYPE_CRNA}=[obj.goodseqs{obj.TYPE_CRNA},randsample(obj.goodseqs{obj.TYPE_T7W},round(goodcleavage*gain*length(obj.goodseqs{obj.TYPE_T7W})),true)];
      obj.active=obj.TYPE_CRNA;
      obj.printdiv(sprintf('T7(gain=%.3f,Good=%.0f%%, Bad=%.0f%%, Bulk=%.0f%%)',gain, goodcleavage*100, badcleavage*100, bulkcleavage*100));
    end

    function RT(obj, efficiency)
    % Reverse transcription
    % bad, good, all get same gain of efficiency
      stopconc=2000;
      obj.cumcost=obj.cumcost+.236*obj.volume;
      rnaconc=obj.conc(obj.TYPE_URNA)+obj.conc(obj.TYPE_CRNA);
      if rnaconc*efficiency>stopconc
        efficiency=stopconc/rnaconc;
        fprintf('Limiting RT efficiency to %.2f due to limited stop oligo\n', efficiency);
      end
      % Apply the gain
      obj.ndna(obj.TYPE_W,:)=obj.ndna(obj.TYPE_W,:)+obj.ndna(obj.TYPE_URNA,:)*efficiency;
      obj.ndna(obj.TYPE_NP,:)=obj.ndna(obj.TYPE_NP,:)+obj.ndna(obj.TYPE_CRNA,:)*efficiency;
      obj.ndna([obj.TYPE_URNA,obj.TYPE_CRNA],:)=0;  % Degrade RNA
      nsamp=round(efficiency*length(obj.goodseqs{obj.TYPE_URNA}));
      if nsamp>0
        obj.goodseqs{obj.TYPE_W}=[obj.goodseqs{obj.TYPE_W},randsample(obj.goodseqs{obj.TYPE_URNA},nsamp,false)];
      end
      nsamp=round(efficiency*length(obj.goodseqs{obj.TYPE_CRNA}));
      if nsamp>0
        obj.goodseqs{obj.TYPE_NP}=[obj.goodseqs{obj.TYPE_NP},randsample(obj.goodseqs{obj.TYPE_CRNA},nsamp,false)];
      end
      obj.goodseqs{obj.TYPE_URNA}=[];
      obj.goodseqs{obj.TYPE_CRNA}=[];
      obj.active=obj.TYPE_NP;
      obj.printdiv(sprintf('RT(eff=%.3f)',efficiency));
    end
    
    function ligate(obj, efficiency)
      nligated=obj.ndna(obj.TYPE_NP,:)*efficiency;
      obj.ndna(obj.TYPE_NP,:)=obj.ndna(obj.TYPE_NP,:)-nligated;
      obj.ndna(obj.TYPE_CIRC,:)=obj.ndna(obj.TYPE_CIRC,:)+nligated;
      samp=randsample(1:length(obj.goodseqs{obj.TYPE_NP}),round(nligated(obj.GOOD)*obj.trackfrac),false);
      nsamp=setdiff(1:length(obj.goodseqs{obj.TYPE_NP}),samp);   % Not circularized ones
      obj.goodseqs{obj.TYPE_CIRC}=[obj.goodseqs{obj.TYPE_CIRC},obj.goodseqs{obj.TYPE_NP}(samp)];
      obj.goodseqs{obj.TYPE_NP}=obj.goodseqs{obj.TYPE_NP}(nsamp);
      obj.active=obj.TYPE_CIRC;
      obj.printdiv(sprintf('Ligate(eff=%.3f)',efficiency));
      obj.cumcost=obj.cumcost+.0002*obj.volume;
    end
    
    function exo(obj,residual)
    % Exo digestion -- remove all DNA but circ
      if nargin<2
        residual=0;
      end
      remove=[obj.TYPE_W,obj.TYPE_NP,obj.TYPE_T7W];
      % Some residual of T7W will appear as W prefix
      obj.ndna(obj.TYPE_W,:)=obj.ndna(obj.TYPE_T7W,:)*residual;
      obj.ndna(obj.TYPE_T7W,:)=0;
      obj.ndna(obj.TYPE_NP,:)=0;
      obj.goodseqs{obj.TYPE_W}=randsample(obj.goodseqs{obj.TYPE_T7W},round(length(obj.goodseqs{obj.TYPE_T7W})*residual),false);
      obj.goodseqs{obj.TYPE_T7W}=0;
      obj.goodseqs{obj.TYPE_NP}=0;
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
      fprintf('%-7.7s%50.50s           %snM                                     ','',note,sprintf('%6.1f ',obs));
      for i=1:length(err)
        rel=err(i);
        str=repmat('*',1,round(rel*10));
        if length(str)>10
          str(10)='|';
        end
        if length(str)>20
          str(20)='|';
        end
        str=['|',str];
        fprintf('%4.2f %-30.30s ',rel, str);
      end
      fprintf('\n');
    end
    
    function qpcr(obj,note,t7,w)
    % Note qPCR measurement (in nM)
      if nargin<4
        w=nan;
      end
      obs=nan(size(obj.ndna,1),1);
      relcirc=sum(obj.ndna(obj.TYPE_CIRC,:))/sum(sum(obj.ndna([obj.TYPE_CIRC,obj.TYPE_T7W],:)));
      obs(obj.TYPE_CIRC)=t7*relcirc;
      obs(obj.TYPE_T7W)=t7*(1-relcirc);
      if isfinite(w)
        obs(obj.TYPE_W)=max(0,w-t7);
      end
      err=t7/nansum(nansum(obj.ndna([obj.TYPE_CIRC,obj.TYPE_T7W],:)));
      if isfinite(w)
        err(2)=w/nansum(nansum(obj.ndna([obj.TYPE_CIRC,obj.TYPE_T7W,obj.TYPE_W],:)));
      end
      err=err*1e-9*obj.volume*1e-6*6.022e23;
      obj.printmeasure(sprintf('qPCR(%.1f,%.1f) %s',t7,w,note),obs,err);
    end
    
    function qubitdna(obj,note,ngul)
    % Qubit HS DNA concentration as given
      ssFactor=1.0;     % Assume ssDNA reads at ssFactor of dsDNA
      factors=[ssFactor ssFactor 1.0 ssFactor 0 0];
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
      factors=[0 0 0 0 1 1];
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
      factors=[1/33 1/33 1/50 1/33 1/40 1/40];
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
      
      % Initial S1 to forward strand (same direction as RNA), S2 as reverse strand
      S1=struct('nbad',sum(obj.ndna([obj.TYPE_T7W],obj.BAD)),...
                'ngood',sum(obj.ndna([obj.TYPE_T7W],obj.GOOD)),...
                'goodseqs',[obj.goodseqs{obj.TYPE_T7W}]);
      S2=struct('nbad',sum(obj.ndna([obj.TYPE_T7W,obj.TYPE_CIRC],obj.BAD)),...
                'ngood',sum(obj.ndna([obj.TYPE_T7W,obj.TYPE_CIRC],obj.GOOD)),...
                'goodseqs',[obj.goodseqs{obj.TYPE_T7W},obj.goodseqs{obj.TYPE_CIRC}]);

      if nargin<3
        maxconc=175;
      end
      maxn=maxconc*1e-9*6.022e23*obj.volume*1e-6;
      for i=1:ncycles
        % Update
        %fprintf('S1=%.2g %.2g S2=%.2g %.2g \n', S1.nbad, S1.ngood, S2.nbad, S2.ngood);
        g1=min(1,(maxn-S1.nbad-S1.ngood)/(S2.nbad+S2.ngood));
        g1=max(g1,0);
        assert(isfinite(g1) && g1>=0);
        S1new=struct('nbad',S1.nbad+S2.nbad*g1,...
                     'ngood',S1.ngood+S2.ngood*g1,...
                     'goodseqs',S1.goodseqs);
        nsamp=round(g1*length(S2.goodseqs));
        if nsamp>0
          S1new.goodseqs=[S1.goodseqs,randsample(S2.goodseqs,nsamp,false)];
        end
        g2=min(1,(maxn-S2.nbad-S2.ngood)/(S1.nbad+S1.ngood));
        g2=max(g2,0);
        assert(isfinite(g2) && g2>=0);
        S2new=struct('nbad',S2.nbad+S1.nbad*g2,...
                     'ngood',S2.ngood+S1.ngood*g2,...
                     'goodseqs',S2.goodseqs);
        nsamp=round(g2*length(S1.goodseqs));
        if nsamp>0
          if length(S1.goodseqs)<=1
            S2new.goodseqs=[S2.goodseqs,S1.goodseqs];
          else
            S2new.goodseqs=[S2.goodseqs,randsample(S1.goodseqs,nsamp,false)];
          end
        end
        S1=S1new;
        S2=S2new;
        %fprintf('Cycle %d: g1=%.2f, g2=%.2f\n', i, g1, g2);
      end
      obj.goodseqs{obj.TYPE_T7W}=[obj.goodseqs{obj.TYPE_T7W},S2.goodseqs];   % Reverse-complement strand is template for future transcription
      %fprintf('S1=%.2g %.2g  S2=%.2g %.2g\n', S1.nbad, S1.ngood, S2.nbad, S2.ngood);
      % No change to number of good sequences since sampling is uniform
      oldtotal=obj.total();
      % During last cycle, assume everything got extended so all parameters are averages
      obj.ndna(obj.TYPE_T7W,obj.GOOD)=obj.ndna(obj.TYPE_T7W,obj.GOOD)+S2.ngood;
      obj.ndna(obj.TYPE_T7W,obj.BAD)=obj.ndna(obj.TYPE_T7W,obj.BAD)+S2.nbad;
      gain=obj.total()/oldtotal;
      obj.cumcost=obj.cumcost+obj.volume*263/250/50;  % Kapa is $263/250U, uses 1U/50ul
      obj.active=obj.TYPE_T7W;
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
