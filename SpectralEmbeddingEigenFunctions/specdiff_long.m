% specdiff


ndata = size(data)

n = prod(ndata)
A = zeros(n,n);
alldata = data(:); % column-wise indexing!

figure; subplot(2,2,1); imagesc(data); title('image','FontSize',16); axis equal;

% pick example points to use throughout script:
disp('pick some points: two close together with very different intensities,');
disp('two neighboring points with very similar intensities, and one non-neighboring');
disp('point with a similar intensity to the previous two.');
[xs(:,1) xs(:,2)] = ginput(5);
xs = round(xs);
xinds = sub2ind(ndata,xs(:,2),xs(:,1));

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
pause;

figure; plot(A(:)); title('absolute distances bw pixels','FontSize',16);
sigma = median(A(:)) %7;% they use 7 for the image of the man in the paper! 
A = exp(-(A).^2/(2*sigma^2));

if 1
     % make A dependent on distance:
    B = zeros(size(data));
    B = reshape(B,ndata);
    [bb(:,1) bb(:,2)] = ind2sub(size(B),[1:n]);
    % check bb:
    figure; scatter(bb(:,1),bb(:,2),[],alldata);
    Adists = zeros(n);
    for j = 1:n
        for k = j:n
            Adists(j,k) = norm(bb(j,:) - bb(k,:),2);
        end
    end
    Adists = Adists + Adists';
    Adists = Adists + diag(ones(n,1)); % make diagonal non-zero
    figure; imagesc(Adists)
    
    A = A./Adists;
    max(A(:)) % check that maximum is along diagonal!!
    figure; imagesc(A)
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
pause;

% 3.2 random walks based on pixel similarity
D = sum(A,2);
D = diag(D);
M = A*inv(D); 

figure; plot(sum(M,1),'m'); title('sum(M,1)','FontSize',20);
figure; plot(ones(1,size(M,1),1)*M,'m'); title('ones*M','FontSize',20);

L = D^(-1/2)*A*D^(-1/2); % = D^(-1/2)*M*D^(1/2)
max(max(L - D^(-1/2)*M*D^(1/2))) % should be near 0!

pause;

% !!!!!!!!!!!
% regularize L, so that it is truly symmetric and all
% eigenvalues/eigenvectors are real:
L = (L + L')/2;
max(max(L - L'))

max(max(L - D^(-1/2)*M*D^(1/2))) % should be near 0!


pause;

[Um,Sm] = eig(M); %eig(L);

[Ul,Sl] = eig(L);
maxs = zeros(n,1);
for k = 1:n
    maxs(k) = max(max(L*Ul(:,k) - Sl(k,k)*Ul(:,k)));
end
figure; plot(maxs)
maxs = zeros(n,1);
Lapprox = D^(-1/2)*M*D^(1/2);
for k = 1:n
    maxs(k) = max(max(Lapprox*Ul(:,k) - Sl(k,k)*Ul(:,k)));
end
figure; plot(maxs)
 
% D^(1/2)*mu are eigenvectors of M with eigenvalues s:
% Lu = su   ... u is a vector, s is a scalar
% D^(-1/2)*M*D^(1/2) u = su
% M*D^(1/2) u = D^(1/2) su
% M*D^(1/2) u = s D^(1/2) u ... x = D^(1/2)*u
% Mx = sx

% sort eigenvectors:
sl = diag(Sl); 
[b ssort] = sort(abs(sl),'descend');
Ul = Ul(:,ssort);
Sl = diag(sl(ssort));
sl = b;

sm = diag(Sm); 
[b ssort] = sort(abs(sm),'descend');
Um = Um(:,ssort);
Sm = diag(sm(ssort));
sm = b;
figure; plot(sm); hold on; plot(sl,'r'); title('eigenvalues of L and M','FontSize',20);


Ulm = D^(1/2)*Ul;
maxs = zeros(n,1);
for k = 1:n
    maxs(k) = max(max(M*Ulm(:,k) - Sl(k,k)*Ulm(:,k)));
end
figure; plot(maxs)

% USE EIGENVECTORS FROM L
S = Sl;
U = Ul;%Ulm;
 
% is Ul orthogonal?
prods = zeros(size(Ul));
for k = 1:size(Ul,2)
    if isreal(Ul(:,k))
        for j =1:size(Ul,2)
            if isreal(Ul(:,j))
                prods(k,j) = Ul(:,k)'*Ul(:,j);
            end
        end
    else
        disp('not real!');
    end
end
figure; imagesc(prods)
figure; plot(diag(prods))
prods(1:10,1:10)
% yay! it's orthogonal
pause;


% make figure 3.3 a -h
nsqu = 7;
figure; 
subplot(nsqu,nsqu,1); imagesc(data); title('image','FontSize',14); axis equal;
for k = 1:nsqu^2-nsqu-1
    subplot(nsqu,nsqu,k+1);
    imagesc(reshape(Um(:,k),ndata)); title(['M eigval ',num2str(k),' = ',num2str(S(k,k))],'FontSize',14);
    axis equal;
end
np = k + 1;
for k = size(U,2) - nsqu + 1:size(U,2) 
    np = np + 1;
    subplot(nsqu,nsqu,np);
    imagesc(reshape(Um(:,k),ndata)); title(['M eigval ',num2str(k),' = ',num2str(S(k,k))],'FontSize',14); 
    axis equal;
end
pause;

% make figure 3.3 a -h with L
nsqu = 7;
figure; 
subplot(nsqu,nsqu,1); imagesc(data); title('image','FontSize',14); axis equal;
for k = 1:nsqu^2-nsqu-1
    subplot(nsqu,nsqu,k+1);
    imagesc(reshape(Ul(:,k),ndata)); title(['L eigval ',num2str(k),' = ',num2str(S(k,k))],'FontSize',14);
    axis equal;
end
np = k + 1;
for k = size(U,2) - nsqu + 1:size(U,2) 
    np = np + 1;
    subplot(nsqu,nsqu,np);
    imagesc(reshape(Ul(:,k),ndata)); title(['L eigval ',num2str(k),' = ',num2str(S(k,k))],'FontSize',14); 
    axis equal;
end



% does the first eigenvector correspond to the stationary distribution?
max(L*Ul(:,1) - Ul(:,1)) % should be near 0
figure; plot(Ul(:,1),'-*'); hold on; plot(L*Ul(:,1),'r'); title('Ul(:,1) (stars) = L*Ul(:,1) (red line)','FontSize',20); 
% yes!!
pause;

% does the first eigenvector correspond to the stationary distribution?
max(M*Ulm(:,1) - Ulm(:,1)) % should be near 0
figure; plot(Ulm(:,1),'-*'); hold on; plot(M*Ulm(:,1),'r'); title('Ulm(:,1) (stars) = M*Ulm(:,1) (red line)','FontSize',20); 


% is U(:,1) positive? yes!
min(U(:,1))

% what about equation 3.5? This need only hold for L, not M:
% "the eigenvector of L associated with the eigenvalue lambda = 1 is given
% by:
alph = sqrt(sum(sum(D)));
d = diag(D);
u1 = (1/alph)*sqrt(d);
figure; plot(u1,'-*'); hold on; plot(L*u1,'r'); % they're identical!
title('u1 (stars) = L*u1 (red line)','FontSize',20); 
figure; plot(u1,'-*'); hold on; plot(Ul(:,1),'r'); % 
title('u1 (stars) = U(:,1) (red line)','FontSize',20); 
max(Ul(:,1) - u1) % should be near 0
% yay!!!!!
% so U(:,1) really is the stationary distribution! (I mean, the eigenvalue
% is 1, so L*u = u,  but I still didn't believe it until now :) )
pause;

figure; subplot(1,3,1); imagesc(reshape(u1,ndata)); title('u1','FontSize',16); axis equal;
subplot(1,3,2); imagesc(reshape(Ulm(:,1),ndata)); title('Ulm(:,1), first eigenvector of M','FontSize',16); axis equal;
subplot(1,3,3); imagesc(reshape(Ul(:,1),ndata)); title('Ul(:,1), first eigenvector of L','FontSize',16); axis equal;
pause;



%% blur kernels and anisotropic diffusion!

% blur = zeros(n,n);
% for j = 1:n % for each pixel
%     ei = zeros(n,1); ei(j) = 1;
%     blur(:,j) = M2*ei;
% end
% blur(:,j) = p_{t,j} = M^{t}*e_{j} = jth column of M_{t}

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


% Reduce dimensions, determine Q to approximate P:

% Choose t and d, so that abs(lambda_{d+1})^t < 1/3...
s = diag(S); % vector of eigenvalues
nt = 100;
st = zeros(size(s,1),nt);
for t = 1:nt
    st(:,t) = abs(s).^t;
end
%hmm....
st(1:30,1:17)
dt = st < 1/3;
dt(1:30,1:15)

d = 14%10%10%15; % take first d eigenvectors
t = 4%8; %4% M^t 1024
Rd = zeros(d,n);
Rd(1:d,1:d) = Rd(1:d,1:d) + diag(ones(d,1));
Ud = U*Rd'; % first d columns of U
isequal(Ud,U(:,1:d))
Sd = S(1:d,1:d); % text says first d columns of the Markov matrix, this uses L...
Dsqrt = D^(1/2);
Dsqrtneg = D^(-1/2);

% These coefficients are used later for the spectral embedding:
wt = zeros(d,n);
pre = ((Sd)^t)*Ud'*Dsqrtneg;
% for x = 1:n
% ei = zeros(n,1); ei(x) = 1;
% wt(:,x) = pre*ei;
% end
wt = pre*diag(ones(size(n)));

Q = zeros(n,n);
pre = Dsqrt*Ud;
% for j = 1:n
%     Q(:,j) = pre*wt(:,j);
% end
Q = pre*wt;
figure; imagesc(Q);

figure;
ind = 3;
x = xinds(ind); subplot(2,2,1); imagesc(reshape(Q(x,:),ndata));  title(['Q ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
hold on; plot(xs(ind,1),xs(ind,2),'go','MarkerSize',10); axis equal;
ind = 4;
x = xinds(4); subplot(2,2,2); imagesc(reshape(Q(x,:),ndata));  title(['Q ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
hold on; plot(xs(ind,1),xs(ind,2),'go','MarkerSize',10); axis equal;
ind = 1;
x = xinds(ind); subplot(2,2,3); imagesc(reshape(Q(x,:),ndata));  title(['Q ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
hold on; plot(xs(ind,1),xs(ind,2),'go','MarkerSize',10); axis equal;
ind = 2;
x = xinds(ind); subplot(2,2,4); imagesc(reshape(Q(x,:),ndata));  title(['Q ',num2str(x), ', t = ',num2str(t)],'FontSize',16);
hold on; plot(xs(ind,1),xs(ind,2),'go','MarkerSize',10); axis equal;

pause;

% get real solution
P = M^t;%M^t;
% Q should approximate P
figure; imagesc(abs(P - Q)); colorbar; title(['Q - P, t = ',num2str(t)],'FontSize',20);
max(max(abs(P-Q)))
pause;
% pretty good!



% try again with higher t!
% t = 32 %... remake Q and P
% % ...
% figure; imagesc(abs(P - Q)); colorbar; title(['Q - P, t = ',num2str(t)],'FontSize',20);
% % lower error, but looks like error is less smooth. Stick with t = 8.
% 
% t = 8;

% 3.4 Spectral Embedding
% each pixel x_i is mapped to the d-dimensional vector wt(i,:).
% wt(i,:).... provides coefficients for the projected blur kernel qt_i
% centered at pixel i. (t = number of times L is diffused (M in text?) )

% p. 49, bottom
% is the first component of wt equal to 1/alpha?
wt(1,:)
1/alph
pause;

Qd = Ud'*D*Ud;
z = (Qd^(1/2))*wt;

% after (3.12): "it follows that q'*q = z'*z
qtq = Q'*Q;
ztz = z'*z;
figure; imagesc(qtq - ztz); colorbar;
max(max(qtq - ztz)) % pretty close!!

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
disp('look at angles between points, and paths/valleys between points');
pause;
% great!

%% ch. 4: SE-MinCut segmentation
n = size(z,2); % number of data points

%% 1) map z to unit sphere
s = zeros(size(z));
for k = 1:n
    s(:,k) = z(:,k)/norm(z(:,k),2);
end
figure; plot(s)
figure; plot(z)
pause;




%% 2) initial guesses for centers

% "inner product is small if pixels belong to different regions, large if
% pixels are near each other in a homogeneous regions"

% start by generating initial guesses for cluster centers:
centers = zeros(n,1); % store indices of centers in here (pre-allocated larger than necessary)
r = floor((n-1)*rand(n,1))+1; % random indices in s
ncenters = 0; % number of centers found so far
nr = 0; % position in r

toofar = [1:n]; % all points that have inner product > tau0 to all centers
% cluster until each point has inner product > tau0 for at least one
% center, so until toofar = []

% initialize first center:
ncenters = ncenters + 1;
nr = nr + 1'
centers(ncenters) = r(nr);
% get inner product of all points to this center:
u = s'*s(:,centers(ncenters));
figure; plot(u);

tau0 = .8%.99;  % .8 in paper.

toofar = find(u < tau0);
size(toofar,1)

while ( (size(toofar,1) > 0) && (nr < n) )
    nr = nr + 1; % look at next point in list of random points
    % "successive samples are constrained to have an inner-product of at
    % least tau0 with regard to previous samples"
    if min(s(:,r(nr))'*s(:,centers(1:ncenters))) < tau0
          % make center
          ncenters = ncenters + 1;
          centers(ncenters) = r(nr);
          % choose clusters until all proints have an inner-product of at
          % least tau0 with some center
          u = s(:,toofar)'*s(:,centers(ncenters));
          toofar = toofar(find(u < tau0));
    end
end
centers = centers(1:ncenters);
ncenters

[a b] = ind2sub(ndata,centers);
figure; imagesc(data);
hold on; plot(b,a,'go','MarkerFaceColor','m','MarkerSize',10);

M = ncenters;
m = s(:,centers);
m'*m

members = zeros(n,M);

%% 3) update cluster centers
maxDiff = 10000;
nrounds = 0;
while maxDiff > .0005
    nrounds = nrounds +1;
    diffs = zeros(M,1);
    for k = 1:M
        oldm = m(:,k);
        wi = max(0,s'*m(:,k) - tau0); % only look at pixels whose inner product with m_k is greater than tau0
        members(find(wi),k) = 1; % record which pixels are used to update this cluster center
        for j = find(wi)'
            m(:,k) = m(:,k) + wi(j)*s(:,j);
        end
        m(:,k) = m(:,k)/sum(wi);
        
        m(:,k) = m(:,k)/norm(m(:,k),2);
        
        diffs(k) = norm(m(:,k) - oldm,2)
        
    end
    maxDiff = max(diffs)
end
nrounds

% Has every pixel been assigned to some cluster?
figure; imagesc(members)
min(max(members,[],2)) % should be 1
% Does every cluster have some members?
max(members)
% how many clusters is each pixel assigned to?
clustersPerpixel = sum(members,2);
figure; plot(clustersPerpixel) % okay...

%% 4) Generate Seed Regions
% S_k := { x_j | s'*m(:,k) >= tau1 }

tau1 = .985;%.999999%.985%.999995; %paper: .985

Seeds = cell(M,1);
nseeds = 6%11; 
for k = 1:M
    u = s'*m(:,k); 
    [a b] = sort(u,'descend');
    Seeds{k} = b(1:nseeds); %find(u >= tau1);%b(1:nseeds); %find(u >= tau1);
 %   a(1:10)
   % pause;
end
Seeds

%% 5) Fig. 4.1: plot seeds
figure; imagesc(data);
Colors = rand(M,3);
hold on;
for k = 1:M
    [a b] = ind2sub(ndata,Seeds{k});
    plot(b,a,'mo','MarkerSize',10,'MarkerFaceColor',Colors(k,:));
   % pause;
end

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

%% 7) eigenfunctions segmentation

% find unique seeds
seeds = cell2mat(Seeds');
uniqseeds = zeros(size(seeds));
nuq = 1;
uniqseeds = seeds;
sunq = size(uniqseeds,2);
k = 1;
while k < sunq
    j = k + 1;
    while j <= sunq
        if isequal(sort(uniqseeds(:,k)),sort(uniqseeds(:,j)))
            uniqseeds = uniqseeds(:,[1:j-1 j+1:end]);
            sunq = sunq - 1;
            j = j - 1;
        end
        j = j + 1;
    end
    k = k + 1;
end
uniqseeds    

addpath('groupmeeting/ssl_eigenfunction');
addpath('groupmeeting');


eigfuncseg(data,uniqseeds(:,[1 2]))
eigfuncseg(data(:,1:42),uniqseeds)


