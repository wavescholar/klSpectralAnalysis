function [Acut bottlenecks] = perturb_easy_im(hh,findxx,sizeim,im)
% does basically the same as itCut, but just trying to understand 
% looking at transitions bw microstates: identifying and cutting bottleneck edges
% following steps in EigenCuts algorithm from Half-lives paper

lag = 2;
A = hh.A;


%% new: try scaling edges by cluster size...
ncl = size(A,1);
clsizes = zeros(ncl,1);
for k = 1:ncl
    clsizes(k) = length(find(hh.cl == k));
end



binNhbr = full(A) > 0;

% initiate
% Step 1:
Ac = A;
ncuts = 100;
niter = 0; % iteration 0
bottlenecks = zeros(prod(size(A)),2);
nnecks = 0;

while ncuts > 0
    niter = niter + 1;
    
    % step 0:
    % take largest connected component of Ac:
    taucc = 0;
    [Acc idividefurther.indices] = cca_wrap(Ac,taucc);
    
    binNhbrcc = binNhbr(idividefurther.indices{1},idividefurther.indices{1});
    clsizescc = clsizes(idividefurther.indices{1});
    clsizesw = clsizescc/sum(clsizescc);
    
    % Step 2:
    dd = sum(Acc,2);
    if max(dd) == 0 % usually means everything's been cut, or was all zeros to begin with
        break; 
    end
    sig = median(dd);
    Dinvsqrt = diag(dd.^(-1/2));
    sqrtD = dd.^(1/2);
    Lc = Dinvsqrt*Acc*Dinvsqrt; % normalized laplacian
    
    % Step 3:
    [U Z V] = svds(Lc,size(Lc,2));
    S = diag(Z);
    
    % plot eigenvectors on ramachandran plot
    if size(Ac,1) == size(U,1) % only easy to do if every cluster is in connected component
        figure(232); 
        for k = 1:min(size(U,2),9)
            subplot(3,3,k);
            uk = U(idividefurther.indices{1}(hh.cl),k);
            ukim = zeros(sizeim);
            ukim(findxx) = uk;
            imagesc(ukim);
            %    scatter(phipsi(1:10:end,1),phipsi(1:10:end,2),[],uk(1:10:end));
            title(['eigval = ',num2str(S(k))]);
        end
        %  pause;
        
    end
    
    
    
    % eigenvalues should decrease in absolute value: abs(s1) > abs(s2) >...
    
    % Step 4:
    % a: Compute implied timescales for each eigenvector
    betak = -lag./log(abs(S));
    % b: compute beta0 in first iteration
    if niter == 1
        epsilon = .25;
        figure; plot(betak(2:end)); title('half-lives of 2nd to final EV','FontSize',20);
        %  x = ginput(1);
        %  beta0 = x(2);
        beta0 = betak(min(length(betak),6));%; % "rough estimate of the smallest half-life one wishes to consider"
        disp(['using beta0 = ',num2str(beta0)]);
    end
    % c: for each eigenvectors with betak > epsilon*beta0, compute the
    % relaxation time sensitivities:
    %   findsens = find(betak > epsilon*beta0);
    findsens = find(betak > beta0);
    
    if isempty(findsens) % end here is no sensitive edges
        break;
    end
    
    lastsens = findsens(end);
    if lastsens == 1 % end here is no sensitive edges
        break;
    end
    
    
    Sij = cell(lastsens,1);
    allsens = [];
    for k = 2:lastsens
        tmp = lag/(S(k) * log(S(k)) * (beta0 * log(S(k)) - lag));
        q = sqrtD.^(-1) .* U(:,k);
        dLogHalf = tmp *(2 * q * q' - ...
            S(k)*q.^2 * ones(1, length(q)) - ...
            S(k)*(ones(length(q),1) *(q.^2)'));
        % fix the diagonal
        Sij{k} = dLogHalf .* binNhbrcc;
        
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % experimental: scale Sij by cluster size:
        %         for j = 1:size(Sij{k},1)
        %            % Sij{k}(j,:) = clsizesw(j)*Sij{k}(j,:);
        %             sSij = sum(Sij{k}(j,:));
        %             if sSij ~= 0
        %                 Sij{k}(j,:) = Sij{k}(j,:)/abs(sSij);
        %             end
        %         end
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        allsens = [allsens ; Sij{k}(find(binNhbrcc))];
    end
    
    maxS = max(max(cell2mat(Sij(2:end))));
    minS = min(min(cell2mat(Sij(2:end))));
    figure(13); clf;
    for k = 2:min(6,lastsens)
        subplot(2,3,k);
        imagesc(Sij{k},[minS maxS]); colorbar;
        title('half-life sensitivities for this eigenvalue','FontSize',16);
    end
    
    
    totBins = 101;
    [ht,bn] = histo(allsens,totBins);
    ht = ht/sum(ht)/abs(bn(1)-bn(2));
    
    figure(231); clf;
    subplot(2,1,1);
    [mm,xx] = max(log10(ht+1));
    grid on;
    plot(bn,log10(ht+1),'k-','linewidth',2);
    hold on; grid on;
    set(gca,'fontsize',20);
    xlabel('dlog(\beta + \beta_0)/d\alpha_{i,j}')
    ylabel('Log Frequency');
    title('all sensitivities of non-zero edges','FontSize',20);
    
    if 0%niter == 1
        %         dLogHalfMaxCut = ginput(1);
        %         dLogHalfMaxCut = dLogHalfMaxCut(1);
        %         disp(['using dLogHalfMaxCut = ',num2str(dLogHalfMaxCut)]);
        x = ginput(1);
        x = x(1);
        cutoff  = x;
        tau = sig*x;
        disp(['using tau = ',num2str(tau)]);
        
    end
    % if 0
    %     [xmax,imax] = extrema(log10(ht+1));
    %     % cutoff should be at a peak smaller than .5*maximum
    %     imaxx = imax(1); % index of maximum
    %     lessxmax = find(xmax < .5*xmax(1));
    %     xmax = xmax(lessxmax);
    %     imax = imax(lessxmax);
    %
    %     lessmax = find(imax < imaxx); % only look at extrema to the left of imax
    %     imax = imax(lessmax);
    %     xmax = xmax(lessmax);
    %
    %     % look at greatest xtrema in lessmax
    %     posscutoffs = bn(imax);
    %
    %     [~,igrlessmax] = sort(posscutoffs,'descend');
    %
    %     if ~isempty(lessmax)
    %         cutoff = bn(imax(igrlessmax(1)));
    %     else
    %         disp('no more lower peaks');
    %         break;
    %     end
    %
    % else
    %     if 1  % alternative: choose first min to the left of the max
    %         % this worked really well, but sometimes chose a very high cutoff
    %         [xmax,imax,xmin,imin] = extrema(log10(ht+1));
    %         imaxx = imax(1); % index of maximum
    %         lessmax = find(imin < imaxx);
    %         xmin = xmin(lessmax);
    %         imin = imin(lessmax);
    %
    % %         % new:
    % %         % only use minima xmin < .5*xmax
    % %         lessxmax = find(xmin < .5*max(xmax));
    % %         xmin = xmin(lessxmax);
    % %         imin = imin(lessxmax);
    %
    %         if ~isempty(imin)
    %             cutoff = bn(max(imin));
    %         else
    %             disp('no more lower peaks');
    %             break;
    %         end
    %     else % try to choose a slightly lower cutoff
    %         [xmax,imax,xmin,imin] = extrema(log10(ht+1));
    %         topmax = find(xmax >= .5*max(xmax))
    %         itopmax = imax(topmax);
    %         imaxx = min(itopmax); % choose maximum that is furthest left
    %         lessmax = find(imin < imaxx);
    %         xmin = xmin(lessmax);
    %         imin = imin(lessmax);
    %
    %         if ~isempty(imin)
    %             cutoff = bn(max(imin));
    %         else
    %             disp('no more lower peaks');
    %             break;
    %         end
    %
    %     end
    %
    % end
    
    if niter == 1
        if 1
            if 1
                [xmax,imax] = extrema(log10(ht+1));
                x = bn(imax(1));
                x = min(x,0); % cutoff shouldn't be larger than 0
                
                
                % cutoff  = x;
                tau = sig*x;%x;%sig*x;
                cutoff  = tau/sig;%x;
                %   tr = cutoff;
                disp(['using tau = ',num2str(tau)]);
            else
                [xmax,imax] = extrema(log10(ht+1));
                % cutoff should be at a peak smaller than .5*maximum
                imaxx = imax(1); % index of maximum
                lessxmax = find(xmax < .5*xmax(1));
                xmax = xmax(lessxmax);
                imax = imax(lessxmax);
                
                lessmax = find(imax < imaxx); % only look at extrema to the left of imax
                imax = imax(lessmax);
                xmax = xmax(lessmax);
                
                % look at greatest xtrema in lessmax
                posscutoffs = bn(imax);
                
                [~,igrlessmax] = sort(posscutoffs,'descend');
                
                if ~isempty(lessmax)
                    cutoff = bn(imax(igrlessmax(1)));
                else
                    disp('no more lower peaks');
                    break;
                end
            end
        else
            if 1  % alternative: choose first min to the left of the max
                % this worked really well, but sometimes chose a very high cutoff
                [xmax,imax,xmin,imin] = extrema(log10(ht+1));
                %   u = find(xmax < .33*max(xmax));
                
                imaxx = imax(1); % index of maximum
                xmaxx = xmax(1); % maximum
                
                % only look at values further left than imaxx
                min_lessmax = find(imin < imaxx);
                max_lessmax = find(imax < imaxx);
                xmin = xmin(min_lessmax);
                imin = imin(min_lessmax);
                xmax = xmax(max_lessmax);
                imax = imax(max_lessmax);
                
                % new:
                % find most-left peak (> .33xmax)
                % only look at minima lefter than this peak (more negative)
                leftpeaks = imax(find(xmax > .5*xmaxx));
                if ~isempty(leftpeaks)
                    imaxx = min(leftpeaks); % furthest left maxima
                else
                    imaxx = imaxx; % use index of maximum
                end
                
                lessxmax = imin(find(imin < imaxx));
                if ~isempty(lessxmax)
                    cutoff = bn(max(lessxmax));
                else
                    disp('no more lower peaks');
                    break;
                end
                
                
                %         % only use minima xmin < .5*xmax
                %         lessxmax = find(xmin < .33*max(xmax));
                %         xmin = xmin(lessxmax);
                %         imin = imin(lessxmax);
                %
                %         if ~isempty(imin)
                %             cutoff = bn(max(imin));
                %         else
                %             disp('no more lower peaks');
                %             break;
                %         end
            else % try to choose a slightly lower cutoff
                [xmax,imax,xmin,imin] = extrema(log10(ht+1));
                topmax = find(xmax >= .5*max(xmax))
                itopmax = imax(topmax);
                imaxx = min(itopmax); % choose maximum that is furthest left
                lessmax = find(imin < imaxx);
                xmin = xmin(lessmax);
                imin = imin(lessmax);
                
                if ~isempty(imin)
                    cutoff = bn(max(imin));
                else
                    disp('no more lower peaks');
                    break;
                end
                
            end
            
            tr = cutoff;%x; % tau/sig, the cutoff used
            
        end
        
    end
    
    %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%%
    % cutoff changes based on sig
    %     cutoff = tau/sig;
    tr = cutoff;
    %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%% %%%%%%%%%%%
    
    ax = axis;
    
    
    plot([tr tr],[ax(3) ax(4)],'k--','linewidth',2)
    
    
    
    % Step 5: non-maximal suppression
    if niter == 1
        binCut0 = zeros(size(Ac,1),size(Ac,1)); % map of edges already cut
    end
    binCut0cc = binCut0(idividefurther.indices{1},idividefurther.indices{1});
    
    binSupp = cell(size(Sij));
    for k = 2:lastsens
     %   binSupp{k} = nonmaxsupp(Sij{k},U(:,k),binNhbrcc,dLogHalfMaxCut,binCut0cc);
     %   binSupp{k} = nonmaxsupp_traj(Sij{k});
     binSupp{k} = zeros(size(binCut0));
        % 1 if suppress, 0 if don't suppress
    end
    
    % Step 6: compute Tk 
 %   tau = -.1;
    T = zeros(1,lastsens);
  %  figure(281); clf;
    for k = 2:lastsens
        Sijtemp = Sij{k};
        
        % compute sum over nonsuppressed edges
        supp = find(binSupp{k});
        Sijtemp(supp) = 0;
        
%         subplot(1,lastsens-1,k-1);
%         [ht,bn] = histo(Sijtemp(:),totBins);
%         ht = ht/sum(ht)/abs(bn(1)-bn(2));
%         [mm,xx] = max(log10(ht+1));
%         plot(bn,log10(ht+1),'k-','linewidth',2);
%         hold on; plot([tr tr],[ax(3) ax(4)],'k--','linewidth',2)
        
        small = find(Sijtemp < cutoff);%tau/sig);
        T(k) = sum(-Sijtemp(small).*Acc(small));
      %  title(['S^{k}ij, T(k) = ',num2str(T(k))],'FontSize',16);
    end
    
    % Step 7: select the eigenmode uk for which Tk is maximal
    % this will be the eigenmode which is response for the most
    % sensitivies*intensity
    [a kstar] = max(T);
    if a == 0
        disp('no more sensitivies < cutoff');
        break;
    end
    
    % Step 8: cut all edges in Ac for which Sij{kstar} < tau/sigma and not
    % suppressed
    % Add intensities back into diagonal.
    % plot:
    % all edges yellow
    figure(231); subplot(2,2,3); 
    spy(A,'g');
    % plot edges remaining before this round in pink:
    subplot(2,2,3); hold on; spy(Ac,'m'); 
    
    Sijtemp = Sij{kstar};
    % remove suppressed edges
    supp = find(binSupp{kstar});
    Sijtemp(supp) = 0;
    % find edges to cut
    small = find(Sijtemp < cutoff);%tau/sig);
    [a b] = ind2sub(size(Acc),small);
    % cut symmetrically:
    smallsym = sub2ind(size(Acc),b,a);
    allsmall = [small; smallsym];
    small = unique(allsmall);
    % ignore indices for which Acc is already 0
    notzero = find(Acc(small) > 0);
    small = small(notzero);
    
    [a b] = ind2sub(size(Acc),small);
    
    ncuts = size(small,1);
    if ncuts > 0
        
        bottlenecks(nnecks+1:nnecks + ncuts,:) = [idividefurther.indices{1}(a) idividefurther.indices{1}(b)];
        nnecks = nnecks + ncuts;
        
        % move density to diagonal:
        precutdd = sum(Acc,2);
        for k = 1:ncuts
            Acc(a(k),a(k)) = Acc(a(k),a(k)) + Acc(a(k),b(k));
        end
        
        cutedges = zeros(ncuts,1);
        for k = 1:ncuts
            cutedges(k) = Acc(a(k),b(k));
        end
        figure(22); subplot(1,2,1); hist(cutedges); title('cut edge values','FontSize',20);
        subplot(1,2,2); hist(Acc(:)); title('all edge values','FontSize',20);
        
        % cut edges:
        Acc(small) = 0;
        disp('is Acc symmetric? (ans should be 0)');
        max(max(Acc - Acc'))
        disp('did degree change? (ans should be 1)');
        postcutdd = sum(Acc,2);%% Set up data for center selection: 

        isequal(precutdd,postcutdd)
        
    end
    
    % embed connected component into Ac
    Ac(idividefurther.indices{1},idividefurther.indices{1}) = Acc;
    
    
    % plot remaining edges in black:
    subplot(2,2,3); hold on; spy(Ac,'k'); axis square;
    title('cut before: green, cut this round: pink','FontSize',18);
    
    % plot edge cuts in rama plot
    % project cuts onto full MD trajectory
    if 0%ncuts > 0
        figure(231); subplot(2,2,4);
        for k = 1:ncuts
            cutsMD = zeros(size(hh.cl,1),1);
            cla = find(hh.cl == a(k));
            clb = find(hh.cl == b(k));
            cutsMD(cla) = 1;
            cutsMD(clb) = 2;
        %    scatter(phipsi(1:10:end,1),phipsi(1:10:end,2),[],cutsMD(1:10:end));
            pause;
        end
    end
    
    % don't re-cut these edges in next round:
    if ncuts > 0
        binCutcc = zeros(size(dd,1),size(dd,1));
        binCutcc(small) = 1;
        binCut = zeros(size(Ac,1),size(Ac,1));
        binCut(idividefurther.indices{1},idividefurther.indices{1}) = binCutcc;
        binCut0 = binCut | binCut0;
    end
  %  pause;
    
end


bottlenecks = bottlenecks(find(bottlenecks(:,1)),:);

Acut = Ac;    
    
 
if 0%nargin > 3 % only include image if hierarchy level 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%  END: Iterative Cutting. %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  Projection Clustering Using K-means %%%%%%%%%%%%%%%%%%%%%
    %%%%  See Ng, Jordan, and Wiess           %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    dd = sum(Ac,2);
    sig = median(dd);
    Dinvsqrt = diag(dd.^(-1/2));
    sqrtD = dd.^(1/2);
    Lc = Dinvsqrt*Ac*Dinvsqrt;
    [U Z V] = svds(Lc,size(Lc,2));
    S = diag(Z);
    betak = -lag./log(abs(S));
    
    minRes2 = 0.25;
    MARKOV = 0; % MARKOV = TRUE; minRes2 = 0.5; nIts = 64;
    
    n = size(Ac,1);
    nIts = niter;
    
    dm = find(abs(betak) < 10000);
    %dm = find(halfL < 1000);
    dm = dm(1) - 1 ;
    
    
    %% Set up data for center selection:
    if MARKOV
        c = U .* (ones(n(1), 1) * (S .^ nIts)');
        E = c;
        E = ((sum(E.*E, 2).^-0.5) * ones(1, size(E,2))) .* E;
    else
        E = U(:, 1:dm);
        E = ((sum(E.*E, 2).^-0.5) * ones(1, dm)) .* E;
    end
    
    %% D0 center selection:
    %% Rerun: rand('state', saveRandState);
    saveRandState = rand('state');
    
    %% Select initial centers to be nearly orthogonal.
    j = round(0.5 + n(1) * rand(1));
    j = max(j, 1); j = min(j,n(1));
    c = E(j,:);
    res = E;
    k = 2;
    while MARKOV | k <= dm
        res = res - (res * (c(k-1,:)')) * c(k-1,:);
        nrmRes2 = sum(res .* res, 2);
        samp = cumsum(nrmRes2 .* (nrmRes2>minRes2));
        if samp(n(1)) == 0 | isnan(samp(n(1))) | isinf(samp(n(1)))
            break;
        end
        samp = samp/samp(n(1));
        r = rand(1);
        idx = find(samp>=r);
        if any(idx)
            c = [c ; E(idx(1), :)];
            k = k+1;
        else
            error('Random draw fanned!??!');
        end
    end
    k = k-1;
    if k < dm & ~MARKOV
        fprintf(2,' Got only %d basis elements\n', k);
        dm = k;
    else
        fprintf(2,' Got %d basis elements\n', k);
        dm = k;
    end
    
    %% Call kmeans
    options = foptions;
    [centers options post errlog] = kmeans(c, E, options);
    discComp = post>0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%% SCRATCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    useImage = 1;
    colourSize = 33;
    fcMap = hsv(colourSize);
    leafComp = discComp;
    figure(503); clf;
    if (useImage)
       % qPts = Pts(:, [2,3]);
       % showIm(reshape(Pts(:,1),sizeIm));
       showIm(im);
       ukim = zeros(sizeim);
            ukim(findxx) = uk;
            imagesc(ukim);
      % qPts = 
    else
      %  qPts = Pts(:, [1,2]);
      %  plot(Pts(:,1), Pts(:,2),'.r');
    end
    hold on;
    for k=1:size(leafComp,2)
        idx = leafComp(:,k);
        if any(idx)
            c = 1+mod(floor(rand*95357),colourSize);
            
            if 0
                plot(qPts(idx,1), qPts(idx>0,2),'o', 'Color', fcMap(c,:),...
                    'MarkerFaceColor', fcMap(c,:));
            end
            ht= text(qPts(idx,1), qPts(idx>0,2),sprintf('%d',k), ...
                'Color', fcMap(c,:), 'FontName', 'Courier', 'FontWeight', 'bold', ...
                'FontSize', 20);
            
            
        else
            fprintf(2,'empty %d component\n', k);
        end
    end
    axis equal; axis off;
    hold off;
    
    
    
    
    
end
    
    
    
    
    
    
    