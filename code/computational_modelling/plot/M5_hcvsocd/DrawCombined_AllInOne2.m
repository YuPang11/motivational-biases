function DrawCombined_AllInOne2()

%% ==============================================================================
%  All-in-one Figure (Width <= 190mm, Right margin = 4mm)  [SVG]
%  Row 1: HC   -> WAIC + All Cues + Go Cues
%  Middle: shared legend (behavior lines + M3c prediction shading)
%  Row 2: OCD  -> WAIC + All Cues + Go Cues
%  Row 3: M5 posteriors (HC vs OCD) 1x5
%
%  Updates kept:
%   (1) Increase spacing between Row2 and Row3
%   (2) Row3: show ylabel + yticks on EVERY subplot
%   (3) Row3 legend: single-line, placed in the blank space between Row2 and Row3
%   (4) Row3 posterior density colors: HC=hc_c_best, OCD=ocd_c_best
%
%  NEW FIXES (this version):
%   (A) Increase overall figure height + increase bottom margin for Row3
%       -> ensure two-line xlabels are fully visible (no clipping)
%
%  FONT + CONSISTENCY FIX (this version):
%   (B) Force ALL interpreters to 'none' (avoid TeX font fallback in SVG)
%   (C) Replace TeX strings with Unicode (Δ, ±, ρ/ε/π/κ) and draw Row3 labels via text()
%   (D) ALL fonts = Arial, size = 9 (truly consistent)
%
%  EXTRA FIX (this version):
%   (E) Move Row3 slightly to the RIGHT to avoid clipping of the left-most subplot
% ==============================================================================

    clc; close all;

    % ======================================================================
    %  (B) FORCE INTERPRETERS TO 'none' (avoid TeX font substitution)
    % ======================================================================
    g = groot;
    oldTI  = get(g,'defaultTextInterpreter');
    oldLTI = get(g,'defaultLegendInterpreter');
    oldATI = get(g,'defaultAxesTickLabelInterpreter');

    try
        set(g,'defaultTextInterpreter','none');
        set(g,'defaultLegendInterpreter','none');
        set(g,'defaultAxesTickLabelInterpreter','none');
    catch
        % older MATLAB may not support some defaults; proceed anyway
    end

    % ensure defaults are restored even if error occurs
    cleanupObj = onCleanup(@() restore_interpreters(oldTI, oldLTI, oldATI));

    % ---------- PATHS ----------
    try
        addpath('F:\PIT\model\helpfunctions'); % boundedline
    catch
    end

    post_csv = 'F:\PIT\modelplot\plot\M5_posteriors_wide.csv';
    save_dir = 'F:\PIT\modelplot\';

    % ---------- GLOBAL LAYOUT (cm) ----------
    figW_cm = 19.0;     % 190mm
    margL   = 1.1;      % left margin
    margR   = 0.8;      % right margin

    % --- Top block (HC/OCD 2x3) ---
    pltH_cm = 4.2;
    margT   = 0.4;
    gapY    = 3.6;      % space between HC and OCD rows (holds shared legend)
    gapX    = 1.20;

    rightScale = 1.12;
    pltW_cm = (figW_cm - margL - margR - 2*gapX) / (2 + rightScale);
    pltW_cm_right = pltW_cm * rightScale;

    % --- Bottom posterior row (Row3) ---
    postH_cm    = 3.8;
    postGapY_cm = 3.0;  % Increase Row2-Row3 gap
    postB_cm    = 2.2;  % increase bottom space
    postGapX_cm = 0.85; % horizontal gaps among the 5 posterior panels

    % ===== Row3 left/right margins =====
    margL_post  = margL + 0.1;
    margR_post  = margR;

    % ===== NEW: move Row3 slightly to the right =====
    postShiftRight_cm = 0.45;          % <-- increase/decrease this if needed
    postRowLeft_cm    = margL_post + postShiftRight_cm;

    % Recompute Row3 width so the whole row still fits inside the figure
    postW_cm = (figW_cm - postRowLeft_cm - margR_post - 4*postGapX_cm) / 5;

    % vertical anchors
    axB_post = postB_cm;
    axB_row2 = axB_post + postH_cm + postGapY_cm;   % OCD row bottom
    axB_row1 = axB_row2 + pltH_cm + gapY;           % HC row bottom

    % increase overall figure height
    figH_cm  = axB_row1 + pltH_cm + margT + 1.2;

    % ---------- TYPOGRAPHY ----------
    fontSz   = 9;
    fontName = 'Arial';
    lwData   = 1.3;
    lwModel  = 2.0;

    % ---------- COLORS ----------
    colWin   = [0 0.60 0.25];
    colAvoid = [0.90 0.20 0.10];

    % HC Theme (Orange)
    hc_c_best  = [229, 157, 66]/255;
    hc_c_other = [246, 228, 198]/255;

    % OCD Theme (Blue)
    ocd_c_best  = [0.25 0.60 0.90];
    ocd_c_other = [0.70 0.85 1.0];

    % ---------- FIGURE ----------
    f = figure('Units','centimeters', ...
        'Position',[5 5 figW_cm figH_cm], ...
        'Color','w', ...
        'PaperPositionMode','auto', ...
        'InvertHardcopy','off');

    set(f,'Renderer','painters'); % SVG friendly

    set(f, 'DefaultAxesFontName', fontName, ...
           'DefaultAxesFontSize', fontSz, ...
           'DefaultTextFontName', fontName, ...
           'DefaultTextFontSize', fontSz);

    set(f,'PaperUnits','centimeters');
    set(f,'PaperSize',[figW_cm figH_cm]);
    set(f,'PaperPosition',[0 0 figW_cm figH_cm]);

    % ---------- GROUP STRUCT ----------
    groups = struct();

    % HC
    groups(1).name     = 'HC';
    groups(1).axBottom = axB_row1;
    groups(1).waic_csv = 'F:\PIT\modelplot\plot\Set0_data.csv';
    groups(1).waic_ref = 'M1';
    groups(1).matfile  = 'F:\PIT\model\EEGmodelHC\HC\results\stan\M5\EEGsingletrialhc3.mat';
    groups(1).csv_pgo  = 'F:\PIT\model\EEGmodelHC\HC\results\stan\M5\M5_osap_samples_pGo.csv';
    groups(1).csv_pco  = 'F:\PIT\model\EEGmodelHC\HC\results\stan\M5\M5_osap_samples_pCorrect.csv';
    groups(1).c_best   = hc_c_best;
    groups(1).c_other  = hc_c_other;

    % OCD
    groups(2).name     = 'OCD';
    groups(2).axBottom = axB_row2;
    groups(2).waic_csv = 'F:\PIT\modelplot\plot\Set0_dataOCD.csv';
    groups(2).waic_ref = 'M1';
    groups(2).matfile  = 'F:\PIT\model\EEGmodelHC\OCD\results\stan\M5\EEGsingletrialocd3.mat';
    groups(2).csv_pgo  = 'F:\PIT\model\EEGmodelHC\OCD\results\stan\M5\M5_osap_samples_pGo.csv';
    groups(2).csv_pco  = 'F:\PIT\model\EEGmodelHC\OCD\results\stan\M5\M5_osap_samples_pCorrect.csv';
    groups(2).c_best   = ocd_c_best;
    groups(2).c_other  = ocd_c_other;

    % ======================================================================
    %  A) TOP BLOCK: HC & OCD (WAIC + Predictions)
    % ======================================================================
    for gidx = 1:2
        G   = groups(gidx);
        axB = G.axBottom;

        % ---------------- Panel 1: WAIC ----------------
        x1  = margL;
        ax1 = axes('Units','centimeters','Position',[x1 axB pltW_cm pltH_cm]); hold(ax1,'on');

        if exist(G.waic_csv,'file')
            T_waic = readtable(G.waic_csv);

            waic_map = containers.Map({'M1','M2','M3','M4','M5'}, {'M1','M2','M3a','M3b','M3c'});
            T_waic.PaperLabel = cell(height(T_waic),1);
            for ii=1:height(T_waic)
                if isKey(waic_map, T_waic.model{ii})
                    T_waic.PaperLabel{ii} = waic_map(T_waic.model{ii});
                else
                    T_waic.PaperLabel{ii} = T_waic.model{ii};
                end
            end

            base_idx = find(strcmp(T_waic.model, G.waic_ref));
            if isempty(base_idx), base_idx = 1; end
            base_val = T_waic.WAIC(base_idx);

            T_waic.Delta  = T_waic.WAIC - base_val;
            T_waic.IsBest = (T_waic.WAIC == min(T_waic.WAIC));

            disp_ord = {'M1','M2','M3a','M3b','M3c'};
            ord_idx = zeros(length(disp_ord),1);
            valid_count = 0;
            for k=1:length(disp_ord)
                ix = find(strcmp(T_waic.PaperLabel, disp_ord{k}));
                if ~isempty(ix)
                    valid_count = valid_count + 1;
                    ord_idx(valid_count) = ix;
                end
            end

            if valid_count > 0
                T_plot = flipud(T_waic(ord_idx(1:valid_count), :));
                nb = height(T_plot);

                for k=1:nb
                    val = T_plot.Delta(k);
                    fc  = T_plot.IsBest(k)*G.c_best + (~T_plot.IsBest(k))*G.c_other;
                    if val <= 0
                        barh(ax1, k, val, 0.65, 'FaceColor', fc, 'EdgeColor','k','LineWidth',0.5);
                    end
                end

                r_row = find(strcmp(T_plot.PaperLabel,'M1'));
                if ~isempty(r_row)
                    line(ax1,[0 0],[r_row-0.35 r_row+0.35],'Color','k','LineWidth',2);
                end

                line(ax1,[0 0],[0.4 nb+0.6],'Color','k','LineWidth',0.5);
                set(ax1,'YLim',[0.4 nb+0.6],'YTick',1:nb,'YTickLabel',T_plot.PaperLabel,'XDir','reverse');

                min_v = min(T_plot.Delta);
                step  = ceil(abs(min_v)/3/50)*50;
                if step<50, step=50; end
                low   = floor(min_v/step)*step;
                xtks  = sort(unique([0:-step:low, 0]));

                set(ax1,'XTick',xtks,'XLim',[min(xtks) 0],'XTickLabelRotation',0);

                xlabel(ax1,'Δ log model evidence (WAIC)', ...
                    'FontSize',fontSz,'FontName',fontName,'Interpreter','none');
            end
        else
            text(ax1,0.5,0.5,'Data Missing','HorizontalAlignment','center', ...
                'FontSize',fontSz,'FontName',fontName,'Interpreter','none');
        end

        set(ax1,'FontSize',fontSz,'FontName',fontName,'Box','off','TickDir','out');

        % ---------------- Load data ----------------
        has_data = exist(G.matfile,'file') && exist(G.csv_pgo,'file') && exist(G.csv_pco,'file');
        if has_data
            S = load(G.matfile);
            data = S.data;
            par  = S.par;

            nSub = par.nSub;
            nRep = 60;
            x_axis = 1:nRep;

            pGo2Win=NaN(nSub,nRep);
            pGo2Avoid=NaN(nSub,nRep);
            pNoGo2Win=NaN(nSub,nRep);
            pNoGo2Avoid=NaN(nSub,nRep);

            pCorrGo2Win=NaN(nSub,nRep);
            pCorrGo2Avoid=NaN(nSub,nRep);

            for iSub=1:nSub
                idx = @(k) data.stim(iSub,:)==k;
                pGo2Win(iSub,:)     = nanmean([data.go(iSub,idx(1)); data.go(iSub,idx(2))],1);
                pGo2Avoid(iSub,:)   = nanmean([data.go(iSub,idx(3)); data.go(iSub,idx(4))],1);
                pNoGo2Win(iSub,:)   = nanmean([data.go(iSub,idx(5)); data.go(iSub,idx(6))],1);
                pNoGo2Avoid(iSub,:) = nanmean([data.go(iSub,idx(7)); data.go(iSub,idx(8))],1);

                pCorrGo2Win(iSub,:)   = nanmean([data.acc(iSub,idx(1)); data.acc(iSub,idx(2))],1);
                pCorrGo2Avoid(iSub,:) = nanmean([data.acc(iSub,idx(3)); data.acc(iSub,idx(4))],1);
            end

            pIncGo2Win   = pGo2Win   - pCorrGo2Win;
            pIncGo2Avoid = pGo2Avoid - pCorrGo2Avoid;

            PGo      = readmatrix(G.csv_pgo);
            PCorrect = readmatrix(G.csv_pco);

            if size(PGo,1)==480, PGo=PGo.'; end
            if size(PCorrect,1)==480, PCorrect=PCorrect.'; end

            pGo_m   = reshape(PGo,      [nSub, 8, nRep]);
            pCorr_m = reshape(PCorrect, [nSub, 8, nRep]);

            Go2Win_m     = squeeze(mean(pGo_m(:,[1 2],:),2));
            Go2Avoid_m   = squeeze(mean(pGo_m(:,[3 4],:),2));
            NoGo2Win_m   = squeeze(mean(pGo_m(:,[5 6],:),2));
            NoGo2Avoid_m = squeeze(mean(pGo_m(:,[7 8],:),2));

            CorrGo2Win_m   = squeeze(mean(pCorr_m(:,[1 2],:),2));
            CorrGo2Avoid_m = squeeze(mean(pCorr_m(:,[3 4],:),2));

            IncGo2Win_m    = Go2Win_m   - CorrGo2Win_m;
            IncGo2Avoid_m  = Go2Avoid_m - CorrGo2Avoid_m;

            bl_plot = @(ax,dat) boundedline(ax, x_axis, mean(dat,1,'omitnan'), ...
                std(dat,0,1,'omitnan')/sqrt(nSub), 'alpha', 'cmap', [0 0 0], 'transparency', 0.65);
        end

        % ---------------- Panel 2: All cues ----------------
        x2  = margL + pltW_cm + gapX;
        ax2 = axes('Units','centimeters','Position',[x2 axB pltW_cm pltH_cm]); hold(ax2,'on');

        if has_data
            plot(ax2,x_axis,mean(pGo2Win,1,'omitnan'),     '-',  'Color',colWin,   'LineWidth',lwData);
            plot(ax2,x_axis,mean(pGo2Avoid,1,'omitnan'),   '-',  'Color',colAvoid, 'LineWidth',lwData);
            plot(ax2,x_axis,mean(pNoGo2Win,1,'omitnan'),   '--', 'Color',colWin,   'LineWidth',lwData);
            plot(ax2,x_axis,mean(pNoGo2Avoid,1,'omitnan'), '--', 'Color',colAvoid, 'LineWidth',lwData);

            [l1,p1]=bl_plot(ax2,Go2Win_m);     set(l1,'Color','k','LineWidth',lwModel); set(p1,'FaceColor',G.c_best,'EdgeColor','none');
            [l2,p2]=bl_plot(ax2,Go2Avoid_m);   set(l2,'Color','k','LineWidth',lwModel); set(p2,'FaceColor',G.c_best,'EdgeColor','none');
            [l3,p3]=bl_plot(ax2,NoGo2Win_m);   set(l3,'Color','k','LineWidth',lwModel); set(p3,'FaceColor',G.c_best,'EdgeColor','none');
            [l4,p4]=bl_plot(ax2,NoGo2Avoid_m); set(l4,'Color','k','LineWidth',lwModel); set(p4,'FaceColor',G.c_best,'EdgeColor','none');
        end

        set(ax2,'YLim',[0 1],'XLim',[1 nRep],'XTick',0:20:60,'Box','off','FontSize',fontSz,'FontName',fontName);
        ylabel(ax2,'P(Go)','FontSize',fontSz,'FontName',fontName,'Interpreter','none');
        xlabel(ax2,'Trial (All cues)','FontSize',fontSz,'FontName',fontName,'Interpreter','none');

        % ---------------- Panel 3: Go cues ----------------
        x3  = margL + 2*pltW_cm + 2*gapX;
        ax3 = axes('Units','centimeters','Position',[x3 axB pltW_cm_right pltH_cm]); hold(ax3,'on');

        if has_data
            plot(ax3,x_axis,mean(pCorrGo2Win,1,'omitnan'), '-',  'Color',colWin,   'LineWidth',lwData);
            plot(ax3,x_axis,mean(pCorrGo2Avoid,1,'omitnan'),'-', 'Color',colAvoid, 'LineWidth',lwData);
            plot(ax3,x_axis,mean(pIncGo2Win,1,'omitnan'),  '--', 'Color',colWin,   'LineWidth',lwData);
            plot(ax3,x_axis,mean(pIncGo2Avoid,1,'omitnan'),'--', 'Color',colAvoid, 'LineWidth',lwData);

            [l1,p1]=bl_plot(ax3,CorrGo2Win_m);   set(l1,'Color','k','LineWidth',lwModel); set(p1,'FaceColor',G.c_best,'EdgeColor','none');
            [l2,p2]=bl_plot(ax3,CorrGo2Avoid_m); set(l2,'Color','k','LineWidth',lwModel); set(p2,'FaceColor',G.c_best,'EdgeColor','none');
            [l3,p3]=bl_plot(ax3,IncGo2Win_m);    set(l3,'Color','k','LineWidth',lwModel); set(p3,'FaceColor',G.c_best,'EdgeColor','none');
            [l4,p4]=bl_plot(ax3,IncGo2Avoid_m);  set(l4,'Color','k','LineWidth',lwModel); set(p4,'FaceColor',G.c_best,'EdgeColor','none');

            if gidx==1
                cX=0.60; cY=0.78;  iX=0.60; iY=0.18;  % HC
            else
                cX=0.58; cY=0.60;  iX=0.58; iY=0.24;  % OCD
            end

            text(ax3, cX, cY, 'Correct Go',   'Units','normalized','FontSize',fontSz,'FontName',fontName,'Interpreter','none');
            text(ax3, iX, iY, 'Incorrect Go', 'Units','normalized','FontSize',fontSz,'FontName',fontName,'Interpreter','none');
        end

        set(ax3,'YLim',[0 1],'XLim',[1 nRep],'XTick',0:20:60,'Box','off','FontSize',fontSz,'FontName',fontName);
        ylabel(ax3,'P(response)','FontSize',fontSz,'FontName',fontName,'Interpreter','none');
        xlabel(ax3,'Trial (Go cues)','FontSize',fontSz,'FontName',fontName,'Interpreter','none');
    end

    % ======================================================================
    %  Shared Legend between HC & OCD rows
    % ======================================================================
    lgd_h = 1.70;
    lgd_y = axB_row2 + pltH_cm + (gapY - lgd_h)/2;
    lgd_left  = margL;
    lgd_width = figW_cm - margL - margR;

    ax_lgd = axes('Units','centimeters', ...
        'Position',[lgd_left, lgd_y, lgd_width, lgd_h], ...
        'XLim',[0 lgd_width], 'YLim',[0 lgd_h], 'Visible','off');
    hold(ax_lgd,'on');

    line_len   = 0.90;
    txt_offset = 0.30;
    row1_y     = lgd_h*0.70;
    row2_y     = lgd_h*0.30;

    col1_x = 0.35;
    col2_x = lgd_width*0.33;
    col3_x = lgd_width*0.60;

    plot(ax_lgd,[col1_x col1_x+line_len],[row1_y row1_y],'-',  'Color',colWin,  'LineWidth',lwData);
    text(ax_lgd,col1_x+line_len+txt_offset,row1_y,'Go-to-Win','FontSize',fontSz,'FontName',fontName,'VerticalAlignment','middle','Interpreter','none');

    plot(ax_lgd,[col2_x col2_x+line_len],[row1_y row1_y],'--','Color',colWin,  'LineWidth',lwData);
    text(ax_lgd,col2_x+line_len+txt_offset,row1_y,'NoGo-to-Win','FontSize',fontSz,'FontName',fontName,'VerticalAlignment','middle','Interpreter','none');

    plot(ax_lgd,[col1_x col1_x+line_len],[row2_y row2_y],'-',  'Color',colAvoid,'LineWidth',lwData);
    text(ax_lgd,col1_x+line_len+txt_offset,row2_y,'Go-to-Avoid','FontSize',fontSz,'FontName',fontName,'VerticalAlignment','middle','Interpreter','none');

    plot(ax_lgd,[col2_x col2_x+line_len],[row2_y row2_y],'--','Color',colAvoid,'LineWidth',lwData);
    text(ax_lgd,col2_x+line_len+txt_offset,row2_y,'NoGo-to-Avoid','FontSize',fontSz,'FontName',fontName,'VerticalAlignment','middle','Interpreter','none');

    patch_h = 0.34;

    patch(ax_lgd, [col3_x col3_x+line_len col3_x+line_len col3_x], ...
        [row1_y-patch_h/2 row1_y-patch_h/2 row1_y+patch_h/2 row1_y+patch_h/2], ...
        groups(1).c_best, 'FaceAlpha',0.65,'EdgeColor','none');
    plot(ax_lgd,[col3_x col3_x+line_len],[row1_y row1_y],'k-','LineWidth',lwModel);
    text(ax_lgd,col3_x+line_len+txt_offset,row1_y,'M3c predictions (mean±sem)', ...
        'FontSize',fontSz,'FontName',fontName,'VerticalAlignment','middle','Interpreter','none');

    patch(ax_lgd, [col3_x col3_x+line_len col3_x+line_len col3_x], ...
        [row2_y-patch_h/2 row2_y-patch_h/2 row2_y+patch_h/2 row2_y+patch_h/2], ...
        groups(2).c_best, 'FaceAlpha',0.65,'EdgeColor','none');
    plot(ax_lgd,[col3_x col3_x+line_len],[row2_y row2_y],'k-','LineWidth',lwModel);
    text(ax_lgd,col3_x+line_len+txt_offset,row2_y,'M3c predictions (mean±sem)', ...
        'FontSize',fontSz,'FontName',fontName,'VerticalAlignment','middle','Interpreter','none');

    % ======================================================================
    %  B) ROW3: M5 posteriors (1x5)
    % ======================================================================
    if exist(post_csv,'file')
        Tp = readtable(post_csv);

        idHC  = strcmp(Tp.Group,'HC');
        idOCD = strcmp(Tp.Group,'OCD');

        post_names = {'rho','epsilon','go_bias','pav_bias','kappa'};
        post_sym_line1 = {'ρ','ε','b','π','κ'};
        post_sym_line2 = {'Feedback sensitivity', 'Learning rate', 'Go bias', ...
                          'Pavlovian bias', 'Instrumental learning bias'};

        col_post_HC  = hc_c_best;
        col_post_OCD = ocd_c_best;
        alpha_post = 0.70;

        y_sym = -0.22;
        y_ann = -0.36;

        for iP = 1:5
            % ===== MODIFIED: Row3 starts from postRowLeft_cm =====
            xP = postRowLeft_cm + (iP-1)*(postW_cm + postGapX_cm);

            axP = axes('Units','centimeters','Position',[xP axB_post postW_cm postH_cm]);
            hold(axP,'on');

            dHC  = Tp.(post_names{iP})(idHC);
            dOCD = Tp.(post_names{iP})(idOCD);

            [maxHC, hdiHC]   = plot_pdf(axP, dHC,  col_post_HC,  alpha_post);
            [maxOCD, hdiOCD] = plot_pdf(axP, dOCD, col_post_OCD, alpha_post);

            base_h = max([maxHC, maxOCD]);
            y_pos  = -0.12*base_h;

            line(axP, hdiHC,  [y_pos y_pos], 'Color', col_post_HC,  'LineWidth', 3);
            line(axP, hdiOCD, [y_pos y_pos], 'Color', col_post_OCD, 'LineWidth', 3);

            axis(axP,'tight');
            yl = ylim(axP);
            ylim(axP, [-yl(2)*0.22, yl(2)*1.08]);

            ylabel(axP,'Density','FontName',fontName,'FontSize',fontSz,'Interpreter','none');
            set(axP,'YTickLabelMode','auto');

            xlabel(axP,'');

            text(axP, 0.5, y_sym, post_sym_line1{iP}, ...
                'Units','normalized','HorizontalAlignment','center', ...
                'FontName',fontName,'FontSize',fontSz,'FontWeight','bold', ...
                'Interpreter','none','Clipping','off');

            text(axP, 0.5, y_ann, post_sym_line2{iP}, ...
                'Units','normalized','HorizontalAlignment','center', ...
                'FontName',fontName,'FontSize',fontSz,'FontWeight','normal', ...
                'Interpreter','none','Clipping','off');

            set(axP,'Box','off','TickDir','out','LineWidth',0.8,'FontSize',fontSz,'FontName',fontName);
            axP.XTickLabelRotation = 0;

            if iP == 2
                xticks(axP,[0.01 0.03]);
                try xtickformat(axP,'%.2f'); catch, end
            elseif iP == 3
                xticks(axP,[0 0.4]);
            elseif iP == 4
                xticks(axP,[-0.5 0.5]);
            end
        end
    else
        warning('Posterior CSV not found: %s (skip posterior row)', post_csv);
    end

    % ======================================================================
    %  Row3 Legend
    % ======================================================================
    lgd3_h = 0.85;
    lgd3_y = axB_post + postH_cm + (postGapY_cm - lgd3_h)/2 - 0.5;
    dy_lgd3_cm = 0.40;
    lgd3_y = lgd3_y - dy_lgd3_cm;

    % ===== MODIFIED: Row3 legend aligns with shifted Row3 =====
    lgd3_left  = postRowLeft_cm;
    lgd3_width = figW_cm - postRowLeft_cm - margR_post;

    ax_lgd3 = axes('Units','centimeters', ...
        'Position',[lgd3_left, lgd3_y, lgd3_width, lgd3_h], ...
        'XLim',[0 lgd3_width], 'YLim',[0 lgd3_h], 'Visible','off');
    hold(ax_lgd3,'on');

    boxW = 0.9;
    boxH = 0.35;
    y0   = lgd3_h/2;

    dx_lgd3_cm = 3;
    xHC = lgd3_width*0.42 + dx_lgd3_cm;
    xO  = lgd3_width*0.58 + dx_lgd3_cm;

    patch(ax_lgd3, [xHC xHC+boxW xHC+boxW xHC], ...
        [y0-boxH/2 y0-boxH/2 y0+boxH/2 y0+boxH/2], ...
        hc_c_best, 'FaceAlpha',0.70, 'EdgeColor',hc_c_best, 'LineWidth',1.0);
    text(ax_lgd3, xHC+boxW+0.25, y0, 'HC', ...
        'FontSize',fontSz,'FontName',fontName,'VerticalAlignment','middle','Interpreter','none');

    patch(ax_lgd3, [xO xO+boxW xO+boxW xO], ...
        [y0-boxH/2 y0-boxH/2 y0+boxH/2 y0+boxH/2], ...
        ocd_c_best, 'FaceAlpha',0.70, 'EdgeColor',ocd_c_best, 'LineWidth',1.0);
    text(ax_lgd3, xO+boxW+0.25, y0, 'OCD', ...
        'FontSize',fontSz,'FontName',fontName,'VerticalAlignment','middle','Interpreter','none');

    % ======================================================================
    %  EXPORT SVG
    % ======================================================================
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end

    out_file = fullfile(save_dir, 'Combined_HC_OCD_plus_M5posteriors_final.svg');

    set(findall(f,'-property','FontName'),'FontName',fontName);
    set(findall(f,'-property','FontSize'),'FontSize',fontSz);

    fprintf('Saving figure to: %s\n', out_file);
    print(f, out_file, '-dsvg', '-r1000');
    fprintf('Done.\n');
end

%% ========================= helpers =========================
function [maxPdf, hdi] = plot_pdf(ax, data, col, alph)
    [pdf_y, pdf_x] = ksdensity(data);
    maxPdf = max(pdf_y);
    hdi = calc_hdi(data, 0.95);

    fill(ax, pdf_x, pdf_y, col, 'FaceAlpha', alph, 'EdgeColor','none');
    plot(ax, pdf_x, pdf_y, 'Color', col, 'LineWidth', 1.2);
end

function hdi_lims = calc_hdi(samples, cred_mass)
    sorted_samples = sort(samples);
    n = length(samples);
    ci_idx = floor(cred_mass * n);
    n_intervals = n - ci_idx;
    interval_widths = sorted_samples(1+ci_idx:n) - sorted_samples(1:n_intervals);
    [~, min_idx] = min(interval_widths);
    hdi_lims = [sorted_samples(min_idx), sorted_samples(min_idx + ci_idx)];
end

function restore_interpreters(oldTI, oldLTI, oldATI)
    g = groot;
    try
        set(g,'defaultTextInterpreter',oldTI);
        set(g,'defaultLegendInterpreter',oldLTI);
        set(g,'defaultAxesTickLabelInterpreter',oldATI);
    catch
    end
end