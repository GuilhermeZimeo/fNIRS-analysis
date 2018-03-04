function BlockAvg(hb, blen, offset, chnlst, keep, plotflag, hbnum)
%Plot block average for each subject (session) provided as input
%
%It is possible to select a subset of channels (chnlst), and/or
%to select the chromophore(s) of interest to be plotted (hbnum)
%
%'keep' is a (#Subj x #Chn) matrix with signal quality labels
%
%The plotflag is a vector (1x3) that allows the user to select
%plots of interest: (1) average per conditions, (2) individual
%responses, (3) all conditions on same graph per chromophore
%
%Note: 'supertitle' function is provided in nirs-toolbox/ext
%
%Last Update: 2018-03-04
%
%By: Guilherme A. Zimeo Morais
%Contact: gazmorais@gmail.com

subjs = length(hb);
conds = length(hb(1).stimulus.keys);
chtot = size(hb(1).data,2)/2;

if nargin < 2
    blen = 5; %baseline prior to onset
end

if nargin < 3
    offset = 10; %display 10s after block ends
end

if nargin < 4 || isempty(chnlst)
    chnlst = 1:chtot;
else
    chtot = length(chnlst);
end

if nargin < 5 || isempty(keep)
    keep = ones(subjs,chtot); %set all channels' quality to 'good'
end

if nargin < 6
    plotflag = [1 0 0]; %only plot average per condition
end

if nargin < 7
    hbnum = [1 2]; %plot both HbO and HbR
end

for n=1:subjs

    Fs = hb(n).Fs;
    wblen = round(Fs*blen);

    hbtyp = unique(hb(n).probe.link.type);
    if length(hbtyp) < length(hbnum)
        hbnum = hbnum(1:length(hbtyp));
    end
    hblen = length(hbnum);
    
    c_dur = zeros(1,conds);
    for c=1:conds
        c_dur(c) = min(hb(n).stimulus.values{c}.dur);
    end
    dur = min(c_dur) + offset;
    
    if plotflag(2)
        fig_t = zeros(1,conds);
        trlim = zeros(conds,2);
    end
    
    if plotflag(3)
        fig_c = zeros(1,hblen);
    end

    for c=1:conds

        trials = length(hb(n).stimulus.values{c}.onset);
        clbl = hb(n).stimulus.keys{c};

        window = length(-wblen:round(dur*Fs));
        block = zeros(trials,window,chtot,hblen);

        twin = (-blen:1/Fs:dur);

        for t=1:trials

            ini = hb(n).stimulus.values{c}.onset(t) - hb(n).time(1);
            fim = ini + dur;

            wmin = round(Fs*ini);
            wmax = round(Fs*fim);  
            wlen = length(wmin-wblen:wmax);

            if wmin > wblen && wmax <= length(hb(n).data)

                for h=1:hblen
                    
                    hbstr = hbtyp{hbnum(h)};
                    hbidx = strcmpi(hb(n).probe.link.type,hbstr);
                    hbdat = hb(n).data(:, hbidx);
                    hbprb = hb(n).probe.link(hbidx, :);

                    base = mean(hbdat( wmin-wblen : wmin , : ),1);

                    for j=1:chtot

                        chn = chnlst(j);
                        chnlbl{j} = [ num2str(hbprb.source(chn)) '-' ...
                                      num2str(hbprb.detector(chn)) ];

                        %check signal quality
                        if keep(n,chn)
                            block(t,(1:wlen),j,h) = hbdat( wmin-wblen : wmax , chn ) - base(chn);
                        else
                            block(t,:,j,h) = NaN;
                        end

                    end

                end

            end

        end
        
        cstr = ['r' 'm'; 'b' 'c'];

        for h = 1:hblen

            curhb = squeeze(block(:,:,:,h));
            hbstr = hbtyp{hbnum(h)};

            if plotflag(1) %plot average per condition
                avg = 1;
                if h==1
                    fig_hb = -1;
                end
                fig_hb = PlotBlockAvg(curhb, twin, cstr(h,c), avg, fig_hb, chnlbl);
                if h==hblen
                    figure(fig_hb);
                    supertitle(['Subj ' num2str(n) '  -  Cond ' clbl '  -  Data All']);
                end
            end

            if plotflag(2) %plot individual responses
                avg = 0;
                fig_t(c) = PlotBlockAvg(curhb, twin, cstr(h,c), avg, -1, chnlbl);
                trlim(c,:) = get(gca, 'YLim');
                figure(fig_t(c));
                supertitle(['Subj ' num2str(n) '  -  Cond ' clbl '  -  Data ' hbstr]);
            end

            if plotflag(3) %plot all conditions on same graph
                avg = 1;
                if c==1
                    fig_c(h) = -1;
                end
                fig_c(h) = PlotBlockAvg(curhb, twin, cstr(h,c), avg, fig_c(h), chnlbl);
                if c==conds
                    figure(fig_c(h));
                    supertitle(['Subj ' num2str(n) '  -  Cond All  -  Data ' hbstr]);
                end
            end

        end                                 

    end

%following is thought to set equal scales for individual plots
%commented out as it makes the supertitle disappear (needs fix)
%     if plotflag(2)
%         trlim = [min(min(trlim)) max(max(trlim))];
%         for h = 1:size(fig_t)
%             fig = figure(fig_t(h));
%             for ax = 1:length(fig.Children)
%                 set(fig.Children(ax), 'YLim', trlim);
%             end
%         end
%     end
    
end

end

function h = PlotBlockAvg(block, twin, cstr, avg, figh, chnlbl)

    if nargin < 3
        cstr = 'r';
    end

    if nargin < 4
        avg = 1;
    end

    if nargin < 5
        figh = -1;
    end

    chtot = size(block,3);

    if figh < 0
        h = figure;
        set(gcf,'Position',[1 41 1920 964],'Color','white');
    else
        h = figh;
        figure(h);
    end

    nRows = floor(sqrt(chtot));
    nCols = ceil(chtot/nRows);

    maxlim = length(twin);

    if avg
        %compute mean
        bmean = squeeze(nanmean(block(:,1:maxlim,:),1));
        %compute standard error of the mean
        bster = squeeze(nanstd(block(:,1:maxlim,:),1)) ./ sqrt(size(block,1));
        %store 'YLim' values
        miny = min(min(bmean-bster)); maxy = max(max(bmean+bster));
    else
        %store all trials
        trials = squeeze(block(:,1:maxlim,:));
        %store 'YLim' values
        miny = min(min(min(trials))); maxy = max(max(max(trials)));
    end

    for ch=1:chtot

        subidx = ch;
        subplot(nRows, nCols, subidx);
        hold on;

        if avg
            shadedErrorBar(twin,squeeze(bmean(:,ch)),squeeze(bster(:,ch)),'lineprops',cstr);
        else
            plot(twin, squeeze(trials(:,:,ch)), cstr);
        end

        if figh < 0 %new figure
            
            set(gca,'XLim',[twin(1) twin(end)]);
            %title(['Channel ' num2str(ch)]);
            title(['Channel ' chnlbl{ch}]);
        
        else   
            
            yref = get(gca,'YLim');
            miny = min(miny, yref(1));
            maxy = max(maxy, yref(2));
            
        end
        
        set(gca,'YLim', [miny maxy]);
        line([0 0],get(gca,'YLim'),'Color','m');

    end

end