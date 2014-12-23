function [M,ZZ]=minconnect(varargin)
%   MINimal CONNECTivity (adjacency) matrix for points on a plane
%   defined by X and Y coordinates (and/or graph of connections). 
%
%   Needs CLINE.
%
%   Applications: graph theory, optimal traffic, astronomy 
%   (e.g. if you want to see a cluster of connected stars 
%   according to certain bounds on distance and/or magnitude etc).
%
%   The connections obey the following optimality condition:
%   breaking any connection divides all points into two groups such that
%   the broken connection corresponds to the shortest distance between
%   the two groups.
%
%   Algorithm: a cluster of already connected points grows by adding 
%   the nearest of resting points
% 
% Call:
%                [M,ZZ]=minconnect(X,Y[,colspec]);% brackets="optional"
%                [M,ZZ]=minconnect(X[,colspec],Y);
%                [M,ZZ]=minconnect(XY[,colspec]);  % XY=[X(:), Y(:)] or XY=X+1i*Y; 
%                 [M,ZZ]=minconnect([colspec,]XY);
%                
%Input:		
%           X = vector of abscissas    
%           Y = vector of ordinates
%           colspec: color/marker/line specification as in CLINE:  
%           if set, connection tree is shown
%   X,Y and colspec (or XY and colspec) may be entered in any
%   sequence, but X should precede Y
%
%Output:	
%      M = minimal connectivity (adjacency) matrix: if points i and j
%       are connected (i<j), then M(i,j)=1 
%       ZZ=(complex) start and finish of all connections;
%DEMO:
%       r=randn(100,2);
%       M = minconnect(r,'2r');figure(gcf),pause
%       MM=minconnect(r(:,1),r(:,2)','bo5');figure(gcf),pause
%       cla;[MMM,ZZ]=minconnect('9g#',r(:,1)+1i*r(:,2));cline(ZZ,'4y');
%       norm(M-MM), norm(M-MMM)

%          Vassili Pastushenko	1998, revised 2006
%==============================================================
R=varargin;
LIN=length(R);
switch LIN>0
    case false
        error('Need coordinates')
    case true  
%Check for plotting
PLIND=false;
   for i=1:LIN
     if ischar(R{i}) 
        PLIND=true;
        colspec=R{i};
        break,
     end
   end
end

if PLIND
    R(i)=[];
end

LIN=length(R);
switch LIN
    case 1
        Z=R{1};
        if size(Z,2)==2
            if ~isreal(Z), 
                error('check coordinates');
            end
        Z=Z(:,1)+1i*Z(:,2);
        else
            Z=Z(:);
            if isreal(Z)|isreal(1i*Z)
                error('Check coordinates')
            end
        end
        
    case 2
        X=R{1};
        Y=R{2};
        Z=X(:)+1i*Y(:);
end

LEN=length(Z);
VER=repmat(Z,1,LEN);
HOR=VER.';
DIST=abs(VER-HOR);
CLUSTER=1;
SOURCE=2:LEN;
M=zeros(LEN);

for i=2:LEN
    WORK=DIST(CLUSTER,SOURCE);
    [mins,ROW]=min(WORK,[],1);
    [dum,COL]=min(mins);
    ROW=ROW(COL);
    INROW=CLUSTER(ROW);
    INSOUR=SOURCE(COL);
    M(INROW,INSOUR)=1;
    CLUSTER(i)=INSOUR;
    SOURCE(COL)=[];
end
M=M+M';
M=M-tril(M);

if PLIND|nargout>1
    [I,J]=find(M);
    x=real(Z);
    y=imag(Z);
    X=[x(I),x(J)];
    Y=[y(I),y(J)];   
    ZZ=(X+1i*Y).';
end

if PLIND
    cline(ZZ,colspec);
    set(gca,'dataaspectratio',[1 1 1])
    figure(gcf)
end
    
