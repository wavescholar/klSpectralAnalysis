function visualize(datasetfname, printflag, outdir)
% visualize(datasetfname, printflag, outdir)
%
% Visualizes the information on Nystrom approximation stored in the file
% datasetfname. If printflag is true, saves the generated graphs to pdfs
% with appropriate names. If outdir is present and printflag is true, the 
% pdfs are saved to this directory.
%
% Graphs generated:
%  - leverage scores and top eigenvalues
%  - spectral, frobenius, and trace norm errors of the exact sampling
%  methods for both the fixed and nonfixed-rank approximants
%  - spectral, frobenius, and trace norm errors of the inexact leverage
%  score sampling methods for both the fixed and nonfixed-rank approximants
%  - timings for the exact sampling methods (only the nonfixed-rank case)
%  - timings for the inexact sampling methods (only the nonfixed-rank case)

%% load the dataset

load(datasetfname);

%% setup plotting

if printflag
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.04 0.03], [0.1 0.04]);
    if ~make_it_tight
        clear subplot;  
    end 

    fontsize = 15;
    width = 6.2;
    height = 12.4;
    basename = fullfile(outdir, savedata.in.datasetbasename);
    printfig = @(figname) printcf([basename figname '.pdf'], fontsize, width, height);
    
    parameterswidth = 6.2;
    parametersheight = 6.2;
    printparameters = @(figname) printcf([basename figname '.pdf'], fontsize, parameterswidth, parametersheight);
    
    timingwidth = 6.2;
    timingheight = 6.2;
    printtiming = @(figname) printcf([basename figname '.pdf'], fontsize, timingwidth, timingheight);
end


%% calculate the mean errors and timings

for lidx = 1:length(savedata.in.lvals)
    % simple
    means = mean(savedata.simpleData(lidx).specerr,2);
    nonfixed_simple_specerr(lidx) = means(1);
    fixed_simple_specerr(lidx) = means(2);
    
    means = mean(savedata.simpleData(lidx).froerr,2);
    nonfixed_simple_froerr(lidx) = means(1);
    fixed_simple_froerr(lidx) = means(2);
    
    means = mean(savedata.simpleData(lidx).trerr,2);
    nonfixed_simple_trerr(lidx) = means(1);
    fixed_simple_trerr(lidx) = means(2);
    
    means = mean(savedata.simpleData(lidx).timings, 2);
    nonfixed_simple_timing(lidx) = means(1);
    fixed_simple_timing(lidx) = means(2);
    
    % srft
    means = mean(savedata.srftData(lidx).specerr,2);
    nonfixed_srft_specerr(lidx) = means(1);
    fixed_srft_specerr(lidx) = means(2);
    
    means = mean(savedata.srftData(lidx).froerr,2);
    nonfixed_srft_froerr(lidx) = means(1);
    fixed_srft_froerr(lidx) = means(2);
    
    means = mean(savedata.srftData(lidx).trerr,2);
    nonfixed_srft_trerr(lidx) = means(1);
    fixed_srft_trerr(lidx) = means(2);

    means = mean(savedata.srftData(lidx).timings, 2);
    nonfixed_srft_timing(lidx) = means(1);
    fixed_srft_timing(lidx) = means(2);
    
    % gaussian
    means = mean(savedata.gaussianData(lidx).specerr,2);
    nonfixed_gaussian_specerr(lidx) = means(1);
    fixed_gaussian_specerr(lidx) = means(2);
    
    means = mean(savedata.gaussianData(lidx).froerr,2);
    nonfixed_gaussian_froerr(lidx) = means(1);
    fixed_gaussian_froerr(lidx) = means(2);
    
    means = mean(savedata.gaussianData(lidx).trerr,2);
    nonfixed_gaussian_trerr(lidx) = means(1);
    fixed_gaussian_trerr(lidx) = means(2);
    
    means = mean(savedata.gaussianData(lidx).timings, 2);
    nonfixed_gaussian_timing(lidx) = means(1);
    fixed_gaussian_timing(lidx) = means(2);
    
    % levscores
    means = mean(savedata.levscoreData(lidx).specerr,2);
    nonfixed_levscore_specerr(lidx) = means(1);
    fixed_levscore_specerr(lidx) = means(2);
    
    means = mean(savedata.levscoreData(lidx).froerr,2);
    nonfixed_levscore_froerr(lidx) = means(1);
    fixed_levscore_froerr(lidx) = means(2);
    
    means = mean(savedata.levscoreData(lidx).trerr,2);
    nonfixed_levscore_trerr(lidx) = means(1);
    fixed_levscore_trerr(lidx) = means(2);
    
    means = mean(savedata.levscoreData(lidx).timings, 2);
    nonfixed_levscore_timing(lidx) = means(1);
    fixed_levscore_timing(lidx) = means(2);
    
    % mixedprobs
%     means = mean(savedata.mixedprobsData(lidx).specerr,2);
%     nonfixed_mixedprobs_specerr(lidx) = means(1);
%     fixed_mixedprobs_specerr(lidx) = means(2);
%     
%     means = mean(savedata.mixedprobsData(lidx).froerr,2);
%     nonfixed_mixedprobs_froerr(lidx) = means(1);
%     fixed_mixedprobs_froerr(lidx) = means(2);
%     
%     means = mean(savedata.mixedprobsData(lidx).trerr,2);
%     nonfixed_mixedprobs_trerr(lidx) = means(1);
%     fixed_mixedprobs_trerr(lidx) = means(2);
%     
%     means = mean(savedata.mixedprobsData(lidx).timings, 2);
%     nonfixed_mixedprobs_timing(lidx) = means(1);
%     fixed_mixedprobs_timing(lidx) = means(2);
    
    % approxlev
    means = mean(savedata.approxlevData(lidx).specerr,2);
    nonfixed_approxlev_specerr(lidx) = means(1);
    fixed_approxlev_specerr(lidx) = means(2);
    
    means = mean(savedata.approxlevData(lidx).froerr,2);
    nonfixed_approxlev_froerr(lidx) = means(1);
    fixed_approxlev_froerr(lidx) = means(2);
    
    means = mean(savedata.approxlevData(lidx).trerr,2);
    nonfixed_approxlev_trerr(lidx) = means(1);
    fixed_approxlev_trerr(lidx) = means(2);
    
    means = mean(savedata.approxlevData(lidx).timings, 2);
    nonfixed_approxlev_timing(lidx) = means(1);
    fixed_approxlev_timing(lidx) = means(2);
    
    % froblev
    means = mean(savedata.froblevData(lidx).specerr,2);
    nonfixed_froblev_specerr(lidx) = means(1);
    fixed_froblev_specerr(lidx) = means(2);
    
    means = mean(savedata.froblevData(lidx).froerr,2);
    nonfixed_froblev_froerr(lidx) = means(1);
    fixed_froblev_froerr(lidx) = means(2);
    
    means = mean(savedata.froblevData(lidx).trerr,2);
    nonfixed_froblev_trerr(lidx) = means(1);
    fixed_froblev_trerr(lidx) = means(2);
    
    means = mean(savedata.froblevData(lidx).timings, 2);
    nonfixed_froblev_timing(lidx) = means(1);
    fixed_froblev_timing(lidx) = means(2);
    
    % speclev
    means = mean(savedata.speclevData(lidx).specerr,2);
    nonfixed_speclev_specerr(lidx) = means(1);
    fixed_speclev_specerr(lidx) = means(2);
    
    means = mean(savedata.speclevData(lidx).froerr,2);
    nonfixed_speclev_froerr(lidx) = means(1);
    fixed_speclev_froerr(lidx) = means(2);
    
    means = mean(savedata.speclevData(lidx).trerr,2);
    nonfixed_speclev_trerr(lidx) = means(1);
    fixed_speclev_trerr(lidx) = means(2);
    
    means = mean(savedata.speclevData(lidx).timings, 2);
    nonfixed_speclev_timing(lidx) = means(1);
    fixed_speclev_timing(lidx) = means(2);
    
end

%% View the leverage scores and spectrums
% for now, not stored to pdf
% lw = 1.5;
% ms = 12;
% maxevaltoshow = max(10, savedata.in.k+1);
% 
% figure();
% stem(savedata.levscores, 'k', 'Linewidth', lw);
% title(sprintf('Leverage scores (k = %d)', savedata.in.k));
% 
% if printflag
%     printparameters('-levscores');
% end
% 
% figure();
% if savedata.in.linearkernelflag
%     s = svds(savedata.in.A, maxevaltoshow).^2;
% else
%     s = svds(savedata.in.A, maxevaltoshow);
% end
% semilogy(s, 'r', 'Linewidth', lw);
% title(sprintf('Top %d eigenvalues', maxevaltoshow));
% 
% if printflag
%     printparameters('-spectrum');
% end

%% Compare the errors of the exact methods
lw = 1.5;
ms = 10;
simple_style = 's-';
simple_color = [.1 .1 .1];
srft_style = 'o-';
srft_color = [.2 .2 .2];
gaussian_style = 'v-';
gaussian_color = [.3 .3 .3];
levscore_style = 'd-';
levscore_color = [.4 .4 .4];
mixedprobs_style = '^-';
mixedprobs_color = [.5 .5 .5];
legendloc = 'Northeast';

figure();
subplot(3,1,1);
plot(savedata.in.lvals, nonfixed_simple_specerr/savedata.optspecerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
hold on;
plot(savedata.in.lvals, nonfixed_srft_specerr/savedata.optspecerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color);
plot(savedata.in.lvals, nonfixed_gaussian_specerr/savedata.optspecerr, ...
     gaussian_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', gaussian_color, 'MarkerFaceColor', gaussian_color);
plot(savedata.in.lvals, nonfixed_levscore_specerr/savedata.optspecerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
%plot(savedata.in.lvals, nonfixed_mixedprobs_specerr/savedata.optspecerr, ...
%     mixedprobs_style, 'LineWidth', lw, 'MarkerSize', ms, ...
%     'MarkerFaceColor', mixedprobs_color);
hold off;
legend('unif', 'srft', 'gaussian', 'levscore', 'Location', legendloc);
title('Relative spectral error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,2);
plot(savedata.in.lvals, nonfixed_simple_froerr/savedata.optfroerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
hold on;
plot(savedata.in.lvals, nonfixed_srft_froerr/savedata.optfroerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color);
plot(savedata.in.lvals, nonfixed_gaussian_froerr/savedata.optfroerr, ...
     gaussian_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', gaussian_color, 'MarkerFaceColor', gaussian_color);
plot(savedata.in.lvals, nonfixed_levscore_froerr/savedata.optfroerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
%plot(savedata.in.lvals, nonfixed_mixedprobs_froerr/savedata.optfroerr, ...
%     mixedprobs_style, 'LineWidth', lw, 'MarkerSize', ms, ...
%     'MarkerFaceColor', mixedprobs_color);
hold off;
legend('unif', 'srft', 'gaussian', 'levscore', 'Location', legendloc);
title('Relative Frobenius error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,3);
plot(savedata.in.lvals, nonfixed_simple_trerr/savedata.opttrerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
hold on;
plot(savedata.in.lvals, nonfixed_srft_trerr/savedata.opttrerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color);
plot(savedata.in.lvals, nonfixed_gaussian_trerr/savedata.opttrerr, ...
     gaussian_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', gaussian_color, 'MarkerFaceColor', gaussian_color);
plot(savedata.in.lvals, nonfixed_levscore_trerr/savedata.opttrerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
%plot(savedata.in.lvals, nonfixed_mixedprobs_trerr/savedata.opttrerr, ...
%     mixedprobs_style, 'LineWidth', lw, 'MarkerSize', ms, ...
%     'MarkerFaceColor', mixedprobs_color);
hold off;
legend('unif', 'srft', 'gaussian', 'levscore', 'Location', legendloc);
title('Relative trace error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

if printflag
    printfig('exact-methods-nonfixed-rank-errors');
end

figure()
subplot(3,1,1);
plot(savedata.in.lvals, fixed_simple_specerr/savedata.optspecerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
hold on;
plot(savedata.in.lvals, fixed_srft_specerr/savedata.optspecerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color);
plot(savedata.in.lvals, fixed_gaussian_specerr/savedata.optspecerr, ...
     gaussian_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', gaussian_color, 'MarkerFaceColor', gaussian_color);
plot(savedata.in.lvals, fixed_levscore_specerr/savedata.optspecerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
%plot(savedata.in.lvals, fixed_mixedprobs_specerr/savedata.optspecerr, ...
%     mixedprobs_style, 'LineWidth', lw, 'MarkerSize', ms, ...
%     'MarkerFaceColor', mixedprobs_color);
hold off;
legend('unif', 'srft', 'gaussian', 'levscore', 'Location', legendloc);
title('Relative spectral error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,2);
plot(savedata.in.lvals, fixed_simple_froerr/savedata.optfroerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
hold on;
plot(savedata.in.lvals, fixed_srft_froerr/savedata.optfroerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color);
plot(savedata.in.lvals, fixed_gaussian_froerr/savedata.optfroerr, ...
     gaussian_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', gaussian_color, 'MarkerFaceColor', gaussian_color);
plot(savedata.in.lvals, fixed_levscore_froerr/savedata.optfroerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
%plot(savedata.in.lvals, fixed_mixedprobs_froerr/savedata.optfroerr, ...
%     mixedprobs_style, 'LineWidth', lw, 'MarkerSize', ms, ...
%     'MarkerFaceColor', mixedprobs_color);
hold off;
legend('unif', 'srft', 'gaussian', 'levscore', 'Location', legendloc);
title('Relative Frobenius error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,3);
plot(savedata.in.lvals, fixed_simple_trerr/savedata.opttrerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
hold on;
plot(savedata.in.lvals, fixed_srft_trerr/savedata.opttrerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color);
plot(savedata.in.lvals, fixed_gaussian_trerr/savedata.opttrerr, ...
     gaussian_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', gaussian_color, 'MarkerFaceColor', gaussian_color);
plot(savedata.in.lvals, fixed_levscore_trerr/savedata.opttrerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
%plot(savedata.in.lvals, fixed_mixedprobs_trerr/savedata.opttrerr, ...
%     mixedprobs_style, 'LineWidth', lw, 'MarkerSize', ms, ...
%     'MarkerFaceColor', mixedprobs_color);
hold off;
legend('unif', 'srft', 'gaussian', 'levscore', 'Location', legendloc);
title('Relative trace error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

if printflag
    printfig('exact-methods-fixed-rank-errors');
end

%% Show the timings for the exact methods
lw = 1.5;
ms = 10;
simple_style = 's-';
simple_color = [.1 .1 .1];
srft_style = 'o-';
srft_color = [.2 .2 .2];
gaussian_style = 'v-';
gaussian_color = [.3 .3 .3];
levscore_style = 'd-';
levscore_color = [.4 .4 .4];
mixedprobs_style = '^-';
mixedprobs_color = [.5 .5 .5];
legendloc = 'Northeast';


figure();
semilogy(savedata.in.lvals, nonfixed_levscore_timing, ...
         levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', levscore_color, 'MarkerFaceColor', levscore_color);
hold on;
semilogy(savedata.in.lvals, nonfixed_simple_timing, ...
         simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', simple_color, 'MarkerFaceColor', simple_color);
semilogy(savedata.in.lvals, nonfixed_srft_timing, ...
         srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', srft_color, 'MarkerFaceColor', srft_color);
semilogy(savedata.in.lvals, nonfixed_gaussian_timing, ...
         gaussian_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', gaussian_color, 'MarkerFaceColor', gaussian_color);
legend('levscore', 'unif', 'srft', 'gaussian', 'Location', legendloc);
xlabel('l');
ylabel('time (s)');

if printflag
    printtiming('exact-methods-timings');
end

%% Compare the errors of the inexact leverage score methods
% display the uniform and leverage sampling errors for calibration
lw = 1.5;
ms = 10;
simple_style = 's-';
simple_color = [.1 .1 .1];
froblev_style = 'o-';
froblev_color = [.2 .2 .2];
approxlev_style = 'v-';
approxlev_color = [.3 .3 .3];
levscore_style = 'd-';
levscore_color = [.4 .4 .4];
speclev_style = '^-';
speclev_color = [.5 .5 .5];
legendloc = 'Northeast';

figure();
subplot(3,1,1);
plot(savedata.in.lvals, nonfixed_levscore_specerr/savedata.optspecerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
hold on;
plot(savedata.in.lvals, nonfixed_simple_specerr/savedata.optspecerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
plot(savedata.in.lvals, nonfixed_approxlev_specerr/savedata.optspecerr, ...
     approxlev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', approxlev_color, 'MarkerFaceColor', approxlev_color);
plot(savedata.in.lvals, nonfixed_froblev_specerr/savedata.optspecerr, ...
     froblev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', froblev_color, 'MarkerFaceColor', froblev_color);
plot(savedata.in.lvals, nonfixed_speclev_specerr/savedata.optspecerr, ...
     speclev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', speclev_color, 'MarkerFaceColor', speclev_color);
hold off;
legend('leverage', 'unif', 'power', 'frob approx lev', ...
       'spectral approx lev', 'Location', legendloc);
title('Relative spectral error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,2);
plot(savedata.in.lvals, nonfixed_levscore_froerr/savedata.optfroerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
hold on;
plot(savedata.in.lvals, nonfixed_simple_froerr/savedata.optfroerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
plot(savedata.in.lvals, nonfixed_approxlev_froerr/savedata.optfroerr, ...
     approxlev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', approxlev_color, 'MarkerFaceColor', approxlev_color);
plot(savedata.in.lvals, nonfixed_froblev_froerr/savedata.optfroerr, ...
     froblev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', froblev_color, 'MarkerFaceColor', froblev_color);
plot(savedata.in.lvals, nonfixed_speclev_froerr/savedata.optfroerr, ...
     speclev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', speclev_color, 'MarkerFaceColor', speclev_color);
hold off;
legend('leverage', 'unif', 'power', 'frob approx lev', ...
       'spectral approx lev', 'Location', legendloc);
title('Relative Frobenius error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,3);
plot(savedata.in.lvals, nonfixed_levscore_trerr/savedata.opttrerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
hold on;
plot(savedata.in.lvals, nonfixed_simple_trerr/savedata.opttrerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
plot(savedata.in.lvals, nonfixed_approxlev_trerr/savedata.opttrerr, ...
     approxlev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', approxlev_color, 'MarkerFaceColor', approxlev_color);
plot(savedata.in.lvals, nonfixed_froblev_trerr/savedata.opttrerr, ...
     froblev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', froblev_color, 'MarkerFaceColor', froblev_color);
plot(savedata.in.lvals, nonfixed_speclev_trerr/savedata.opttrerr, ...
     speclev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', speclev_color, 'MarkerFaceColor', speclev_color);
hold off;
legend('leverage', 'unif', 'power', 'frob approx lev', ...
       'spectral approx lev', 'Location', legendloc);
title('Relative trace error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

if printflag
    printfig('inexact-methods-nonfixed-rank-errors');
end

figure();
subplot(3,1,1);
plot(savedata.in.lvals, fixed_levscore_specerr/savedata.optspecerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
hold on;
plot(savedata.in.lvals, fixed_simple_specerr/savedata.optspecerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
plot(savedata.in.lvals, fixed_approxlev_specerr/savedata.optspecerr, ...
     approxlev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', approxlev_color, 'MarkerFaceColor', approxlev_color);
plot(savedata.in.lvals, fixed_froblev_specerr/savedata.optspecerr, ...
     froblev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', froblev_color, 'MarkerFaceColor', froblev_color);
plot(savedata.in.lvals, fixed_speclev_specerr/savedata.optspecerr, ...
     speclev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', speclev_color, 'MarkerFaceColor', speclev_color);
hold off;
legend('leverage', 'unif', 'power', 'frob lev', ...
       'spectral lev', 'Location', legendloc);
title('Relative spectral error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,2);
plot(savedata.in.lvals, fixed_levscore_froerr/savedata.optfroerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
hold on;
plot(savedata.in.lvals, fixed_simple_froerr/savedata.optfroerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
plot(savedata.in.lvals, fixed_approxlev_froerr/savedata.optfroerr, ...
     approxlev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', approxlev_color, 'MarkerFaceColor', approxlev_color);
plot(savedata.in.lvals, fixed_froblev_froerr/savedata.optfroerr, ...
     froblev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', froblev_color, 'MarkerFaceColor', froblev_color);
plot(savedata.in.lvals, fixed_speclev_froerr/savedata.optfroerr, ...
     speclev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', speclev_color, 'MarkerFaceColor', speclev_color);
hold off;
legend('leverage', 'unif', 'power', 'frob lev', 'spectral lev', 'Location', legendloc);
title('Relative Frobenius error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,3);
plot(savedata.in.lvals, fixed_levscore_trerr/savedata.opttrerr, ...
     levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', levscore_color, 'MarkerFaceColor', levscore_color);
hold on;
plot(savedata.in.lvals, fixed_simple_trerr/savedata.opttrerr, ...
     simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', simple_color, 'MarkerFaceColor', simple_color);
plot(savedata.in.lvals, fixed_approxlev_trerr/savedata.opttrerr, ...
     approxlev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', approxlev_color, 'MarkerFaceColor', approxlev_color);
plot(savedata.in.lvals, fixed_froblev_trerr/savedata.opttrerr, ...
     froblev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', froblev_color, 'MarkerFaceColor', froblev_color);
plot(savedata.in.lvals, fixed_speclev_trerr/savedata.opttrerr, ...
     speclev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', speclev_color, 'MarkerFaceColor', speclev_color);
hold off;
legend('leverage', 'unif', 'power', 'frob lev', ...
       'spectral lev', 'Location', legendloc);
title('Relative trace error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

if printflag
    printfig('inexact-methods-fixed-rank-errors');
end

%% Show the timings for the inexact methods
lw = 1.5;
ms = 10;
simple_style = 's-';
simple_color = [.1 .1 .1];
froblev_style = 'o-';
froblev_color = [.2 .2 .2];
approxlev_style = 'v-';
approxlev_color = [.3 .3 .3];
levscore_style = 'd-';
levscore_color = [.4 .4 .4];
speclev_style = '^-';
speclev_color = [.5 .5 .5];
legendloc = 'Northeast';

figure();
semilogy(savedata.in.lvals, nonfixed_levscore_timing, ...
         levscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', levscore_color, 'MarkerFaceColor', levscore_color);
hold on;
semilogy(savedata.in.lvals, nonfixed_simple_timing, ...
         simple_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', simple_color, 'MarkerFaceColor', simple_color);
semilogy(savedata.in.lvals, nonfixed_approxlev_timing, ...
         approxlev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', approxlev_color, 'MarkerFaceColor', approxlev_color);
semilogy(savedata.in.lvals, nonfixed_froblev_timing, ...
         froblev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', froblev_color, 'MarkerFaceColor', froblev_color);
semilogy(savedata.in.lvals, nonfixed_speclev_timing, ...
         speclev_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', speclev_color, 'MarkerFaceColor', speclev_color);
legend('levscore', 'unif', 'power', 'frob lev', ...
       'spectral lev', 'Location', legendloc);
xlabel('l');
ylabel('time (s)');

if printflag
    printtiming('inexact-methods-timings');
end

%% close all figures
if printflag
  %  close all
end

end
