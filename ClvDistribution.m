classdef ClvDistribution < matlab.mixin.Copyable
  properties
    seqs;   % Vector of sequences showing sampling
    cleavage;   % Cleavage of each sequence id
  end
  
  methods
    function obj=ClvDistribution(cleavage,nseqs)
    % Build the distribution using nseqs ids with given cleavage distribution
      assert(nseqs<=length(cleavage));
      obj.seqs=1:nseqs;
      obj.cleavage=cleavage;
    end

    function clear(obj)
      obj.seqs=[];
    end
    
    function n=nseq(obj)
      n=length(obj.seqs);
    end
    
    function [nu,n]=nunique(obj)
    % Compute the mean number of copies of each good sequence
      nu=length(unique(obj.seqs));
      n=length(obj.seqs);
    end

    function c=meanCleavage(obj)
      c=mean(obj.cleavage(obj.seqs));
    end
    
    function resample(obj,n,replacement)
      if n==0
        obj.seqs=[];
      else
        obj.seqs=randsample(obj.seqs,n,replacement);
      end
    end
    
    function f=frachigh(obj,mincleavage)
    % Return fraction of pool that has mincleavage or higher
      clv=obj.cleavage(unique(obj.seqs));
      f=mean(clv>=mincleavage);
    end

    function k=khigh(obj,mincleavage)
    % Return avg number of copies of high cleavers
      hseqs=obj.seqs(obj.cleavage(obj.seqs)>=mincleavage);
      k=length(hseqs)/length(unique(hseqs));
    end

    function addsample(obj,src,n,replacement)
      if n>0
        obj.seqs=[obj.seqs,randsample(src.seqs,n,replacement)];
      end
    end
    
    function [u,c]=splitBasedOnCleavage(obj)
    % Split into 2 new distributions based on random cleavage
      clv=obj.cleavage(obj.seqs);
      selclv=rand(size(clv))<clv;
      u=ClvDistribution(obj.cleavage,0);
      c=ClvDistribution(obj.cleavage,0);
      u.seqs=obj.seqs(~selclv);
      u.cleavage=obj.cleavage;
      c.seqs=obj.seqs(selclv);
      c.cleavage=obj.cleavage;
    end

    function [a,b]=splitRandom(obj,p)
    % Split into 2 new distributions based on random cleavage
      sela=rand(size(obj.seqs))<p;
      a=ClvDistribution(obj.cleavage,0);
      b=ClvDistribution(obj.cleavage,0);
      a.seqs=obj.seqs(sela);
      a.cleavage=obj.cleavage;
      b.seqs=obj.seqs(~sela);
      b.cleavage=obj.cleavage;
    end
    
    function plotdist(obj,ti)
    % Plot distribution of counts of sequences
      if nargin<2
        ti='Divtrack.dist';
      end
      cnt=hist(obj.seqs,1:max(obj.seqs));
      [dist,n]=hist(cnt,0:max(cnt));
      dist=dist/sum(dist)*obj.cdist{active}.nunique();
      setfig(ti);clf;
      plot(n(n>0),dist(n>0),'o-');
      xlabel('Number of copies of sequence');
      ylabel('Number of sequences');
      title(sprintf('%s: Total of %.2f sequences with mean of %.2f copies',ti, obj.cdist{active}.nunique(), obj.avgcopies()));
    end

  end
end

