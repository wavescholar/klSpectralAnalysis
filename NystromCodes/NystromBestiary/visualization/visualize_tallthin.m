function visualize_tallthin(datasetfname1, datasetfname2, printflag, outdir)
% visualize(datasetfname1, datasetfname2, printflag, outdir)
%
% Visualizes the information on Nystrom approximation stored in the files
% datasetfname1 and datasetfname2. If printflag is true, saves the generated graphs to pdfs
% with appropriate names. If outdir is present and printflag is true, the 
% pdfs are saved to this directory.
%
% datasetfname1 is assumed to contain the information on all inexact
% leverage score methods other than the one using the tall thin
% approximation of leverage scores
%
% datasetfname2 is assumed to contain only the information on the tall
% thin approximation algorithm
%
% Graphs generated:
%  - spectral, frobenius, and trace norm errors of the inexact leverage
%  score sampling methods for both the fixed and nonfixed-rank approximants
%  - timings for the inexact sampling methods (only the nonfixed-rank case)

%% setup plotting
load(datasetfname1) % so basename can be defined

if printflag
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.05], [0.045 0.03], [0.1 0.04]);
    if ~make_it_tight
        clear subplot;  
    end 

    fontsize = 15;
    width = 6.2;
    height = 12.4;
    basename = fullfile(outdir, savedata.in.datasetbasename);
    printfig = @(figname) printcf([basename figname '.pdf'], fontsize, width, height);
    
    timingwidth = 6.2;
    timingheight = 6.2;
    printtiming = @(figname) printcf([basename figname '.pdf'], fontsize, timingwidth, timingheight);
end


%% calculate the mean errors and timings

load(datasetfname1);

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
    
    % leverage
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
end

load(datasetfname2);

for lidx = 1:length(savedata.in.lvals)
    % tall thin lev score approx alg
    means = mean(savedata.tallthinData(lidx).specerr,2);
    nonfixed_tallthin_specerr(lidx) = means(1);
    fixed_tallthin_specerr(lidx) = means(2);
    
    means = mean(savedata.tallthinData(lidx).froerr,2);
    nonfixed_tallthin_froerr(lidx) = means(1);
    fixed_tallthin_froerr(lidx) = means(2);
    
    means = mean(savedata.tallthinData(lidx).trerr,2);
    nonfixed_tallthin_trerr(lidx) = means(1);
    fixed_tallthin_trerr(lidx) = means(2);
    
    means = mean(savedata.tallthinData(lidx).timings, 2);
    nonfixed_tallthin_timing(lidx) = means(1);
    fixed_tallthin_timing(lidx) = means(2);
    
    % QR-computed levscores
    means = mean(savedata.qrlevscoreData(lidx).specerr, 2);
    nonfixed_qrlevscore_specerr(lidx) = means(1);
    fixed_qrlevscore_specerr(lidx) = means(2);
    
    means = mean(savedata.qrlevscoreData(lidx).froerr, 2);
    nonfixed_qrlevscore_froerr(lidx) = means(1);
    fixed_qrlevscore_froerr(lidx) = means(2);
    
    means = mean(savedata.qrlevscoreData(lidx).trerr, 2);
    nonfixed_qrlevscore_trerr(lidx) = means(1);
    fixed_qrlevscore_trerr(lidx) = means(2);
    
    means = mean(savedata.qrlevscoreData(lidx).timings, 2);
    nonfixed_qrlevscore_timing(lidx) = means(1);
    fixed_qrlevscore_timing(lidx) = means(2);
    
end


%% Compare the errors of the inexact leverage score methods
% display the uniform and QR leverage sampling errors for calibration
lw = 1.5;
ms = 10;
simple_style = 's-';
simple_color = [.1 .1 .1];
froblev_style = 'o-';
froblev_color = [.2 .2 .2];
approxlev_style = 'v-';
approxlev_color = [.3 .3 .3];
qrlevscore_style = 'd-';
qrlevscore_color = [.4 .4 .4];
speclev_style = '^-';
speclev_color = [.5 .5 .5];
tallthin_style = 'x-';
tallthin_color = [.6 .6 .6];
srft_style = '+-'; 
srft_color = [.7 .7 .7];

legendloc = 'Northeast';

figure();
subplot(3,1,1);
plot(savedata.in.lvals, nonfixed_qrlevscore_specerr/savedata.optspecerr, ...
     qrlevscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', qrlevscore_color, 'MarkerFaceColor', qrlevscore_color);
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
plot(savedata.in.lvals, nonfixed_srft_specerr/savedata.optspecerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color, 'Color', srft_color);
plot(savedata.in.lvals, nonfixed_tallthin_specerr/savedata.optspecerr, ...
     tallthin_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', tallthin_color, 'MarkerFaceColor', tallthin_color);
hold off;
legend('QR lev', 'unif', 'power', 'frob approx lev', ...
       'spectral approx lev', 'srft', 'alg 1', 'Location', legendloc);
title('Relative spectral error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,2);
plot(savedata.in.lvals, nonfixed_qrlevscore_froerr/savedata.optfroerr, ...
     qrlevscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', qrlevscore_color, 'MarkerFaceColor', qrlevscore_color);
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
plot(savedata.in.lvals, nonfixed_srft_froerr/savedata.optfroerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color, 'Color', srft_color);
plot(savedata.in.lvals, nonfixed_tallthin_froerr/savedata.optfroerr, ...
     tallthin_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', tallthin_color, 'MarkerFaceColor', tallthin_color);
hold off;
legend('QR lev', 'unif', 'power', 'frob approx lev', ...
       'spectral approx lev', 'srft', 'alg 1', 'Location', legendloc);
title('Relative Frobenius error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,3);
plot(savedata.in.lvals, nonfixed_qrlevscore_trerr/savedata.opttrerr, ...
     qrlevscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', qrlevscore_color, 'MarkerFaceColor', qrlevscore_color);
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
plot(savedata.in.lvals, nonfixed_srft_trerr/savedata.opttrerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color, 'Color', srft_color);
plot(savedata.in.lvals, nonfixed_tallthin_trerr/savedata.opttrerr, ...
     tallthin_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', tallthin_color, 'MarkerFaceColor', tallthin_color);
hold off;
legend('QR lev', 'unif', 'power', 'frob approx lev', ...
       'spectral approx lev', 'srft', 'alg 1', 'Location', legendloc);
title('Relative trace error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

if printflag
    printfig('inexact-methods-plus-tallthin-nonfixed-rank-errors');
end

figure();
subplot(3,1,1);
plot(savedata.in.lvals, fixed_qrlevscore_specerr/savedata.optspecerr, ...
     qrlevscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', qrlevscore_color, 'MarkerFaceColor', qrlevscore_color);
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
plot(savedata.in.lvals, fixed_srft_specerr/savedata.optspecerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color, 'Color', srft_color);
plot(savedata.in.lvals, fixed_tallthin_specerr/savedata.optspecerr, ...
     tallthin_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', tallthin_color, 'MarkerFaceColor', tallthin_color);
hold off;
legend('QR lev', 'unif', 'power', 'frob lev', ...
       'spectral lev', 'srft', 'alg 1', 'Location', legendloc);
title('Relative spectral error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,2);
plot(savedata.in.lvals, fixed_qrlevscore_froerr/savedata.optfroerr, ...
     qrlevscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', qrlevscore_color, 'MarkerFaceColor', qrlevscore_color);
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
plot(savedata.in.lvals, fixed_srft_froerr/savedata.optfroerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color, 'Color', srft_color);
plot(savedata.in.lvals, fixed_tallthin_froerr/savedata.optfroerr, ...
     tallthin_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', tallthin_color, 'MarkerFaceColor', tallthin_color);
hold off;
legend('QR lev', 'unif', 'power', 'frob lev', ...
       'spectral lev', 'srft', 'alg 1', 'Location', legendloc);
title('Relative Frobenius error');
extents = axis();
xlabel('l');
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

subplot(3,1,3);
plot(savedata.in.lvals, fixed_qrlevscore_trerr/savedata.opttrerr, ...
     qrlevscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', qrlevscore_color, 'MarkerFaceColor', qrlevscore_color);
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
plot(savedata.in.lvals, fixed_srft_trerr/savedata.opttrerr, ...
     srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', srft_color, 'MarkerFaceColor', srft_color, 'Color', srft_color);
plot(savedata.in.lvals, fixed_tallthin_trerr/savedata.opttrerr, ...
     tallthin_style, 'LineWidth', lw, 'MarkerSize', ms, ...
     'Color', tallthin_color, 'MarkerFaceColor', tallthin_color);
hold off;
legend('QR lev', 'unif', 'power', 'frob lev', ...
       'spectral lev', 'srft', 'alg 1', 'Location', legendloc);
title('Relative trace error');
xlabel('l');
extents = axis();
axis([extents(1) extents(2) max(extents(3), 0), extents(4)])

if printflag
    printfig('inexact-methods-plus-tallthin-fixed-rank-errors');
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
qrlevscore_style = 'd-';
qrlevscore_color = [.4 .4 .4];
speclev_style = '^-';
speclev_color = [.5 .5 .5];
tallthin_style = 'x-';
tallthin_color = [.6 .6 .6];
srft_style = '+-'; 
srft_color = [.7 .7 .7];

legendloc = 'Southeast';

figure();
semilogy(savedata.in.lvals, nonfixed_qrlevscore_timing, ...
         qrlevscore_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', qrlevscore_color, 'MarkerFaceColor', qrlevscore_color);
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
semilogy(savedata.in.lvals, nonfixed_srft_timing, ...
         srft_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', srft_color, 'MarkerFaceColor', srft_color, 'Color', srft_color);
semilogy(savedata.in.lvals, nonfixed_tallthin_timing, ...
         tallthin_style, 'LineWidth', lw, 'MarkerSize', ms, ...
         'Color', tallthin_color, 'MarkerFaceColor', tallthin_color);
legend('QR lev', 'unif', 'power', 'frob lev', ...
       'spectral lev', 'srft', 'alg 1', 'Location', legendloc);
xlabel('l');
ylabel('time (s)');

if printflag
    printtiming('inexact-methods-plus-tallthin-timings');
end

%% close all figures
if printflag
    close all
end

end