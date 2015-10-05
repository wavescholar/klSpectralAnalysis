% specdiff_mem.m 
% trying to be more memory/time efficient to work with large images


%% get information about data
ndata = size(data)

n = prod(ndata)
alldata = data(:); % column-wise indexing!

figure; subplot(2,2,1); imagesc(data); title('image','FontSize',16); axis equal;

% pick example points to use throughout script:
disp('pick some points: two close together with very different intensities,');
disp('two neighboring points with very similar intensities, and one non-neighboring');
disp('point with a similar intensity to the previous two.');
[xs(:,1) xs(:,2)] = ginput(5);
xs = round(xs);
xinds = sub2ind(ndata,xs(:,2),xs(:,1));

if 0
    % A is spatially dependent:
   dists = squareform(pdist(alldata).^2); % square bc pdist gives norm, and I want regular distance
   sigma = 10%median(dists(:));
   window = 5;
   A = brightAffty(data,3,1);
end
    

A = squareform(pdist(alldata).^2); % square bc pdist gives norm, and I want regular distance

subplot(2,2,2); imagesc(A); title('absolute distances bw pixels','FontSize',16); axis equal;
% A makes more sense if you look at the distances from individual pixels:
ind = 1;
subplot(2,2,3); imagesc(reshape(A(:,xinds(ind)),ndata)); title('distance from one pixel','FontSize',16);
hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerFaceColor','m','MarkerSize',10); axis equal;
subplot(2,2,1); hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerFaceColor','m','MarkerSize',10); axis equal;
ind = 3;
subplot(2,2,4); imagesc(reshape(A(:,xinds(ind)),ndata)); title('distance from another pixel','FontSize',16);
hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerFaceColor','m','MarkerSize',10); axis equal;
subplot(2,2,1); hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerFaceColor','m','MarkerSize',10); axis equal;
pause(5);

%% make affinity matrix
sigma = median(A(:)) %7;% they use 7 for the image of the man in the paper! 
A = exp(-(A).^2/(2*sigma^2));

if 1
    disp('there must be a way to make B with a look-up table!!');
    disp('also I can sparsify to make both A and B more efficiently');
     % make A dependent on distance:
    B = zeros(size(data));
    B = reshape(B,ndata);
    [bb(:,1) bb(:,2)] = ind2sub(size(B),[1:n]);
    clear B;
    % check bb:
  %  figure; scatter(bb(:,1),bb(:,2),[],alldata);
    Adists = zeros(n);
    for j = 1:n
        if mod(j,100) == 1
            j
        end
        for k = j:n
            Adists(j,k) = norm(bb(j,:) - bb(k,:),2);
        end
    end
    Adists = Adists + Adists';
    Adists = Adists + diag(ones(n,1)); % make diagonal non-zero
%    figure; imagesc(Adists)
    
    A = A./Adists;
    clear Adists;
    max(A(:)) % check that maximum is along diagonal!!
 %   figure; imagesc(A)
end

figure; subplot(2,2,1); imagesc(data); title('image','FontSize',16); axis equal;
subplot(2,2,2); imagesc(A); title('affinities between pixels','FontSize',16); axis equal;
% A makes more sense if you look at the distances from individual pixels:
ind = 1;%471;%1 %90; 
subplot(2,2,3); imagesc(reshape(A(:,xinds(ind)),ndata)); title('affinites from outer border pixel','FontSize',16);
hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerFaceColor','m','MarkerSize',10); axis equal;
subplot(2,2,1); hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerFaceColor','m','MarkerSize',10); axis equal;
ind = 3; %471; 
subplot(2,2,4); imagesc(reshape(A(:,xinds(ind)),ndata)); title('affinites from inner square pixel','FontSize',16);
hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerFaceColor','m','MarkerSize',10); axis equal;
subplot(2,2,1); hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerFaceColor','m','MarkerSize',10); axis equal;
pause(5);


%% 3.2 random walks based on pixel similarity
dd = sum(A,2);
% D = diag(dd);
invD = diag(dd.^(-1));
% never actually use M: (just L)
%M = A*invD;  

%figure; plot(sum(M,1),'m'); title('sum(M,1)','FontSize',20);
%figure; plot(ones(1,size(M,1),1)*M,'m'); title('ones*M','FontSize',20);

Dinvsqrt = diag(dd.^(-1/2));
Dsqrt = diag(dd.^(1/2));
L = Dinvsqrt*A*Dinvsqrt;
%L = D^(-1/2)*A*D^(-1/2); % = D^(-1/2)*M*D^(1/2)
% max(max(L - Dinvsqrt*M*Dsqrt)) % should be near 0!

% !!!!!!!!!!!
% regularize L, so that it is truly symmetric and all
% eigenvalues/eigenvectors are real:
L = (L + L')/2;
max(max(L - L'))

pause(1);


%% USE EIGENVECTORS FROM L
% try out lanczos method next!!!
if n > 2500
    tic
    disp('calculating top 50 eigenvalues:');
    [U,S] = jdqr(L,min(50,n));
    toc
else
    tic
    disp('calculating all eigenvalues:');
    [U,S] = eig(L);
    toc
end

maxs = zeros(size(U,2),1);
for k = 1:size(U,2)
    maxs(k) = max(max(L*U(:,k) - S(k,k)*U(:,k)));
end
figure; plot(maxs); title('how true are the eigenvectors: LU = sU?'); 

% sort eigenvectors:
s = diag(S); 
[b ssort] = sort(abs(s),'descend');
U = U(:,ssort);
S = diag(s(ssort));
s = s(ssort); % not = b, because b is absolute values..
figure; plot(s,'m'); title('eigenvalues of L and M','FontSize',20);


% make figure 3.3 a -h with L
nsqu = 7;
figure; 
subplot(nsqu,nsqu,1); imagesc(data); title('image','FontSize',14); axis equal;
for k = 1:nsqu^2-nsqu-1
    subplot(nsqu,nsqu,k+1);
    imagesc(reshape(U(:,k),ndata)); title(['L eigval ',num2str(k),' = ',num2str(S(k,k))],'FontSize',14);
    axis equal;
end
np = k + 1;
for k = size(U,2) - nsqu + 1:size(U,2) 
    np = np + 1;
    subplot(nsqu,nsqu,np);
    imagesc(reshape(U(:,k),ndata)); title(['L eigval ',num2str(k),' = ',num2str(S(k,k))],'FontSize',14); 
    axis equal;
end

% is U(:,1) positive? yes!
disp(['is the first eigenvector positive? min = ',num2str(min(U(:,1)))]);
if min(U(:,1)) < 0
    U(:,1) = -U(:,1);
    min(U(:,1))
end


% does the first eigenvector correspond to the stationary distribution?
max(L*U(:,1) - U(:,1)) % should be near 0
figure(301); subplot(3,1,1); plot(U(:,1),'-*'); hold on; plot(L*U(:,1),'r'); title('U(:,1) (stars) = L*U(:,1) (red line)','FontSize',16); 
% yes!!
% % does the first eigenvector correspond to the stationary distribution?
% max(M*Ulm(:,1) - Ulm(:,1)) % should be near 0
% figure; plot(Ulm(:,1),'-*'); hold on; plot(M*Ulm(:,1),'r'); title('Ulm(:,1) (stars) = M*Ulm(:,1) (red line)','FontSize',20); 

% what about equation 3.5? This need only hold for L, not M:
% "the eigenvector of L associated with the eigenvalue lambda = 1 is given
% by:
alph = sqrt(sum(dd));
u1 = (1/alph)*sqrt(dd);
figure(301); subplot(3,1,2); plot(u1,'-*'); hold on; plot(L*u1,'r'); % they're identical!
title('u1 (stars) = L*u1 (red line)','FontSize',16); 
figure(301); subplot(3,1,3); plot(u1,'-*'); hold on; plot(U(:,1),'r'); % 
title('u1 (stars) = U(:,1) (red line)','FontSize',16); 
max(U(:,1) - u1) % should be near 0
% yay!!!!!
% so U(:,1) really is the stationary distribution! (I mean, the eigenvalue
% is 1, so L*u = u,  but I still didn't believe it until now :) )
pause(2);

figure; subplot(1,3,1); imagesc(reshape(u1,ndata)); title('u1','FontSize',16); axis equal;
% subplot(1,3,2); imagesc(reshape(Ulm(:,1),ndata)); title('Ulm(:,1), first eigenvector of M derived from L','FontSize',16); axis equal;
subplot(1,3,3); imagesc(reshape(U(:,1),ndata)); title('U(:,1), first eigenvector of L','FontSize',16); axis equal;
pause(2);



%% blur kernels and anisotropic diffusion!

disp(' not viewing blur kernels because time-consuming. use specdiff_long.m to view these!');
if 0
    Mdiff = M;%M;
    t = 1;
    nsmooth = 25;
    nsm = floor(sqrt(nsmooth));
    
    for k = 1:nsmooth
        t = t*2;
        Mdiff = Mdiff*Mdiff;
        P = Mdiff;
        
        figure(223);
        ind = 3;
        x = xinds(ind); subplot(2,2,1); imagesc(reshape(P(x,:),ndata));  title(['P ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
        hold on; plot(xs(ind,1),xs(ind,2),'go','MarkerSize',10); axis equal;
        ind = 4; %471;
        x = xinds(ind); subplot(2,2,2); imagesc(reshape(P(x,:),ndata));  title(['P ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
        hold on; plot(xs(ind,1),xs(ind,2),'go','MarkerSize',10); axis equal;
        ind = 1;
        x = xinds(ind); subplot(2,2,3); imagesc(reshape(P(x,:),ndata));  title(['P ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
        hold on; plot(xs(ind,1),xs(ind,2),'go','MarkerSize',10); axis equal;
        ind = 2;
        x = xinds(ind); subplot(2,2,4); imagesc(reshape(P(x,:),ndata));  title(['P ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
        hold on; plot(xs(ind,1),xs(ind,2),'go','MarkerSize',10); axis equal;
        
        % anisotropically smoothed image:
        B = zeros(1,n);
        for j = 1:n
            B(j) = alldata'*P(:,j);
        end
        maxx = max([max(B) max(alldata)]);
        minn = min([min(B) min(alldata)]);
        figure(224); subplot(1,2,1); imagesc(reshape(data,ndata)); title('image','FontSize',12); axis equal;
        subplot(1,2,2); imagesc(reshape(B,ndata)); title(['anisotropically diffused image, t = ',num2str(t)],'FontSize',16); axis equal;
        
        % look at inner products between blur kernels:
        % inner product small if pixels belong to different regions, and large if
        % pixels belong to the same region
        innblur = zeros(n,n);
        for kk = 1:n
            for j = 1:n
                innblur(kk,j) = P(:,kk)'*P(:,j);
            end
        end
        figure(11); subplot(nsm,nsm,k); imagesc(innblur); title('inner products bw blur kernels - gets smaller for larger t','FontSize',12); colorbar;
        max(innblur(:)) - min(innblur(:))
        
        % Note: as t increases, all blur kernels approach the stationary
        % distribution
        f = mean(P(:,1)./u1)
        figure(14); subplot(nsm,nsm,k); plot(P); hold on; plot(f*u1,'r') % maybe...
        title('as t increases, all blur kernels approach the stat dist','FontSize',12);
        
        pause;
    end
end

%% Reduce dimensions, determine Q to approximate P:

% Choose t and d, so that abs(lambda_{d+1})^t < 1/3...
nt = 30; % number of times to diffuse M... could make smaller, since I only ever look at the top 30
st = zeros(size(s,1),nt);
abss = abs(s);
for t = 1:nt
    st(:,t) = abss.^t;
end
%hmm....
disp('actually want to choose t so that M^t has not yet reached its stationary');
disp('distribution... otherwise all the z vectors will be basically the same.');
disp('text says to choose s.y. lambda_{d+1}^t < 1/3, so look for largest d,t combo with 1/3 as value');
st(1:min(size(st,1),30),1:end)
% dt = st < 1/3;
% dt(1:30,1:15)

d = input(' number of eigenvectors, d = ?');
t = input(' number of diffusions, t = ?');

if 0 % used this for sparse A
    d = min(d,size(st,1))
    t = 21
end

% Rd = zeros(d,n);
% Rd(1:d,1:d) = Rd(1:d,1:d) + diag(ones(d,1));
Ud = U(:,1:d);
Sd = S(1:d,1:d); % text says first d columns of the Markov matrix, this uses L...
sd = s(1:d);

% These coefficients are used later for the spectral embedding:
Sdt = diag(sd.^t); % this is sigItsDiag in getSeeds.m
pre = Sdt*Ud';
wt = pre*Dinvsqrt; clear pre; % wt is Ek in getSeeds.m
%clear Dinvqrt;

% p. 49, bottom
% is the first component of wt equal to 1/alpha?
figure; plot(wt(1,:),'m');
hold on; plot(round(n/2),1/alph,'kx','MarkerSize',10); 
title('is the first component of wt (pink line) = 1/alpha (x)','FontSize',20);
pause(2);
disp('if this is wrong, check that S is correct!') 


if 0 % never actually need this, but I like the figs
    % projected blur kernels
    pre = Dsqrt*Ud;
    %clear Dsqrt;
    Q = pre*wt;
    if size(Q,1) < 5000
        figure; imagesc(Q);
    end
    
    figure;
    ind = 3;
    x = xinds(ind); subplot(2,2,1); imagesc(reshape(Q(x,:),ndata));  title(['Q ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
    hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerSize',10); axis equal;
    ind = 4;
    x = xinds(4); subplot(2,2,2); imagesc(reshape(Q(x,:),ndata));  title(['Q ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
    hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerSize',10); axis equal;
    ind = 1;
    x = xinds(ind); subplot(2,2,3); imagesc(reshape(Q(x,:),ndata));  title(['Q ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
    hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerSize',10); axis equal;
    ind = 2;
    x = xinds(ind); subplot(2,2,4); imagesc(reshape(Q(x,:),ndata));  title(['Q ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
    hold on; plot(xs(ind,1),xs(ind,2),'mo','MarkerSize',10); axis equal;
end

pause(5);


%% 3.4 Spectral Embedding
% each pixel x_i is mapped to the d-dimensional vector wt(i,:).
% wt(i,:).... provides coefficients for the projected blur kernel qt_i
% centered at pixel i. (t = number of times L is diffused (M in text?) )


%Qd = Ud'*D*Ud;
pre = Ud'*diag(dd); 
Qd = pre*Ud;  % this is Q in getSeeds.m
% see how chakra does this - more lines, but maybe more efficient?
z = (Qd^(1/2))*wt; % this is Zk in getSeeds.m

% % after (3.12): "it follows that q'*q = z'*z
% qtq = Q'*Q;
% ztz = z'*z;
% figure; imagesc(qtq - ztz); colorbar;
% max(max(qtq - ztz)) % pretty close!!

% look at fig. 3.7, c (but with different combinations of components of z, instead of just z_2 and z_4)
figure; 
n = 0;
x = xinds(1)% 206; % look at two neighboring points with VERY different intensities
y = xinds(2)%207; 
f1 = xinds(3); % look at two neighboring points with similar intensities
f2 = xinds(4); 
f3 = xinds(5); % and one non-neighboring point with a similar intensity to f1 and f2
nd = round(sqrt(min(7,d)*(min(7,d)/2)))
for k = 1:min(7,d)
    for j = k:min(7,d)
        n = n + 1;
        subplot(nd,nd+1,n); 
        scatter(z(k,:),z(j,:),'yo');
        xlabel(['z_',num2str(k)]);
        ylabel(['z_',num2str(j)]);
        hold on; plot(z(k,x),z(j,x),'mx');
        hold on; plot(z(k,y),z(j,y),'gx');
        hold on; plot(z(k,f1),z(j,f1),'bx');
        hold on; plot(z(k,f2),z(j,f2),'bx');
        hold on; plot(z(k,f3),z(j,f3),'kx');
    end
    set(gca,'FontSize',12);
end
subplot(nd,nd+1,3); title('red and green should be separated','FontSize',16);
disp('look at angles between points, and paths/valleys between points');
disp('If these are all linear, choose a different t,d combo - lower t, higher d?');
pause(5);
% great!

%% ch. 4: SE-MinCut segmentation

%[s m seeds uniqseeds] = zseg(z,data);
[s m seeds uniqseeds] = getSeeds_eigfunc(z,data);
M = size(seeds,2);


%% 6) Naive Segmentation
u = s'*m; 
[minm lbl] = min(u,[],2);

figure; imagesc(data);
Colors = rand(M,3);
hold on;
for k = 1:M
    cl = find(lbl == k);
    [a b] = ind2sub(ndata,cl);
    plot(b,a,'ms','MarkerSize',10,'MarkerFaceColor',Colors(k,:));
  %  pause;
end
title('Naive Segmentation','FontSize',20);

%% 7) eigenfunctions segmentation

% Try to make data square for now, because eigfuncs isn't very adaptable
% right now. 
if size(data,2) ~= size(data,1)
    DATA = data(1:min(ndata),1:min(ndata));
    se = zeros(size(uniqseeds));
    for k = 1:size(se,2)
        [a b] = ind2sub(ndata,uniqseeds(:,k));
        for j = 1:size(se,1)
            if a(j) < size(DATA,1)
                if b(j) < size(DATA,2)
                    se(j,k) = sub2ind(size(DATA),a(j),b(j));
                end
            end
        end
    end
else
    DATA = data;
    se = uniqseeds;
end
    
NUM_EVECS = 15; % use 15 eigenfunctions
figs = 0; % don't plot figures
[dd2, uu2] = eigfuncs(DATA,NUM_EVECS,figs);
f_efunc = funcseg(DATA,se,dd2,uu2);

%% ADD NEW SEEDS! 
if 1 
    % add new seeds
    M = size(se,2);
    figure; imagesc(data);
    Colors = rand(M+1,3);
    hold on;
    for k = 1:M
        [a b] = ind2sub(ndata,se(:,k));
        plot(b,a,'mo','MarkerSize',10,'MarkerFaceColor',Colors(k,:));
    end
    newseeds = floor(ginput(size(se,1)));
    newseedssub = (newseeds(:,1)-1)*size(data,1) + newseeds(:,2);
    % add points to plot:
    k = k + 1;
    [a b] = ind2sub(ndata,newseedssub);
    plot(b,a,'mo','MarkerSize',10,'MarkerFaceColor',Colors(k,:));
    
    se = [se newseedssub];
    
    %% redo segmentation with new seeds!
    f_efunc = funcseg(DATA,se,dd2,uu2);
end
    













