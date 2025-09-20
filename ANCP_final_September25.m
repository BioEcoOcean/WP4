% ANCP_By_Float_Per_Year_with_Diagnostics.m
%
% PURPOSE
%   Estimate Annual Net Community Production (ANCP) from BGC-Argo nitrate profiles
%   by detecting a winter maximum and late-summer minimum in depth-integrated nitrate,
%   converting the drawdown (N) to carbon (C) with the Redfield ratio, and producing
%   diagnostics (time series, maps, histograms, and profile evolution plots).
%
% KEY NOTES
%   • Vertical domain: 0–150 m on a 1 m grid (common_depth).
%   • Interpolation is linear in depth; pressure is made strictly monotonic via unique().
%   • Seasonal windows: winter max = Dec(y?1) & Jan–Mar(y); summer min = Aug–Nov(y).
%   • Units: nitrate (µmol/kg) is converted to mol/m^3
%     using ? (kg/m^3) and then integrated over meters; result is in mol/m^2 and finally
%     converted to gC/m^2 via molar mass of carbon (12 g/mol) and Redfield C:N (106/16).
%   • Selected floats (selected_floats) will get detailed per-float time-series plots.
%
% REQUIREMENTS
%   MATLAB with Mapping Toolbox for geoscatter/geobasemap (optional but used here).
%   A folder 'argo_data_nitrate' with one Argo profile per NetCDF file.

R = 106 / 16;                          % Redfield C:N ratio (mol C per mol N)
common_depth = 0:1:150;                % Target regular depth grid [m]
depth_idx = find(common_depth <= 150); % Indices for the integration range (explicit)
ref_time = datetime(1950,1,1);         % CF/Argo epoch for JULD
data_dir = 'argo_data_nitrate';        % Directory with *.nc files (one profile per file)
nc_files = dir(fullfile(data_dir, '*.nc')); % List all NetCDF inputs

% Parse unique float IDs from filenames (e.g., 'BD6903549_165.nc' -> 'BD6903549')
float_names = unique(cellfun(@(x) regexp(x,'^[A-Z]+\d{7}','match','once'),{nc_files.name},'UniformOutput',false));
float_names = float_names(~cellfun('isempty',float_names));

% Aggregated outputs across floats/years
all_lat = []; all_lon = []; all_ANCP = []; all_years = [];
ANCP_profile_lats = {}; ANCP_profile_lons = {}; all_integrated_chla = []; float_ids = {};
chlorophyll_valid_counts = containers.Map; % Count valid chlorophyll integrals per float

% Store time–depth profile stacks and integrals for diagnostics (per float)
float_evolution = struct('float_id', {}, 'dates', {}, 'nitrate_profiles', {}, 'chlorophyll_profiles', {}, 'integrated_nitrate', {}, 'integrated_chla', {});

% Floats for which to produce detailed per-float time series (left/right axes)
selected_floats = {'BD6901485','BD5903887','SD6903550','SD6903574','SD6901485','SD5903887'};

% ===============================
% MAIN LOOP: process each float ID
% ===============================
for f = 1:length(float_names)
    float_id = float_names{f};
    float_files = dir(fullfile(data_dir,[float_id,'*.nc'])); % All files for this float
    profiles = struct([]);                                   % Interpolated profiles live here
    fprintf('Processing float: %s\n',float_id);

    % ---------------------------------------------
    % Read and parse each NetCDF file for the float
    % ---------------------------------------------
    for k = 1:length(float_files)
        try
            file_path = fullfile(data_dir,float_files(k).name);
            info = ncinfo(file_path);
            vars = {info.Variables.Name};

            % Resolve variable names (supporting common aliases)
            nitrate_var  = intersect(vars,{'NITRATE_ADJUSTED','NTAW_ADJUSTED'});
            pressure_var = intersect(vars,{'PRES','PRESSURE'});
            chla_var     = intersect(vars,{'CHLA_FLUORESCENCE_ADJUSTED','CHLA_ADJUSTED'});
            time_var     = intersect(vars,{'JULD','TIME'});
            lat_var      = intersect(vars,{'LATITUDE'});
            lon_var      = intersect(vars,{'LONGITUDE'});

            % Require presence of core variables to proceed
            if isempty(nitrate_var)||isempty(pressure_var)||isempty(time_var)||isempty(lat_var)||isempty(lon_var)
                continue;
            end

            % Read arrays
            nitrate   = ncread(file_path,nitrate_var{1});    % µmol/kg, [depth × profileColumn]
            pressure  = ncread(file_path,pressure_var{1});   % dbar (? m), [depth × profileColumn]
            juld      = ncread(file_path,time_var{1});       % days since 1950-01-01
            time      = ref_time + days(juld);               % convert to datetime
            latitude  = ncread(file_path,lat_var{1});        % deg N
            longitude = ncread(file_path,lon_var{1});        % deg E
            chla = [];
            if ~isempty(chla_var), chla = ncread(file_path,chla_var{1}); end % mg/m^3 or similar

            % ---------------------------------------------------------
            % Iterate over nitrate profile columns (each is a cast/time)
            % ---------------------------------------------------------
            for i = 1:size(nitrate,2)
                % Column-wise pressure/nitrate and validity mask (0–150 m)
                Pn_col = pressure(:,i); N_col = nitrate(:,i);
                idxN   = ~isnan(Pn_col)&~isnan(N_col)&Pn_col<=150;
                if sum(idxN)<2, continue; end                      % Need ?2 valid points to interpolate
                Pn = Pn_col(idxN); N = N_col(idxN);
                if numel(unique(Pn))<2, continue; end              % interp1 needs monotonic (unique) x

                % Depth-interpolate nitrate to common grid (enforce monotonic via unique)
                try
                    [Puniq, iu] = unique(Pn);                      % strictly increasing pressure
                    Nuniq = N(iu);
                    if numel(Puniq) >= 2
                        N_interp = interp1(Puniq,Nuniq,common_depth,'linear',NaN);
                    else
                        continue;
                    end
                catch, continue; end

                % Prepare chlorophyll interpolation; keep the best-covered column (most non-NaNs)
                CH_interp = NaN(size(common_depth));
                if ~isempty(chla) && ~all(all(isnan(chla)))
                    best_count=-Inf;                                % sentinel for coverage comparison
                    for ch=1:size(chla,2)
                        CH_col = chla(:,ch);
                        % Use the corresponding pressure column if it exists; otherwise skip
                        if ch<=size(pressure,2), Pc_col=pressure(:,ch); else, continue; end
                        % Only consider CHL columns sampled on the same grid length as pressure
                        if ~isequal(size(Pc_col),size(CH_col)), continue; end
                        % Valid paired P–CHL samples in 0–150 m
                        CH_idx = ~isnan(Pc_col)&~isnan(CH_col)&Pc_col<=150;
                        if sum(CH_idx)>=2
                            CH_p=Pc_col(CH_idx); CH_v=CH_col(CH_idx);
                            try
                                [CHuniq, iu_ch] = unique(CH_p);    % enforce monotonic x
                                CHuniq_v = CH_v(iu_ch);
                                if numel(CHuniq) >= 2
                                    ch_try = interp1(CHuniq,CHuniq_v,common_depth,'linear',NaN);
                                    % Keep the interpolation with the most valid (non-NaN) depth points
                                    cur_count=sum(~isnan(ch_try(depth_idx)));
                                    if cur_count>best_count
                                        best_count=cur_count; CH_interp=ch_try;
                                    end
                                end
                            catch, continue; end
                        end
                    end
                end

                % Store this cast's interpolated fields and metadata
                profiles(end+1).nitrate=N_interp;
                profiles(end).chlorophyll=CH_interp;
                profiles(end).date=time(min(i,end));        % time stamp for this column
                profiles(end).latitude=latitude(1);         % per-file lat/lon metadata
                profiles(end).longitude=longitude(1);
            end
        catch ME
            % Report and continue on file-level errors (keeps pipeline robust)
            fprintf('Error reading %s: %s\n',float_files(k).name,ME.message);
        end
    end

    % Skip floats that produced no usable profiles
    if isempty(profiles), continue; end

    % ------------------------------------------------------
    % Assemble time series and compute depth-integrated fields
    % ------------------------------------------------------
    nitrate_interp_all=NaN(length(common_depth),length(profiles));       % [depth × time]
    chlorophyll_interp_all=NaN(length(common_depth),length(profiles));    % [depth × time]
    integrated_nitrate=NaN(1,length(profiles));                           % ? nitrate dz (0–150 m)
    integrated_chla=NaN(1,length(profiles));                              % ? chlorophyll dz (0–150 m)
    dates_all=NaT(1,length(profiles));                                    % time stamps
    lats=NaN(1,length(profiles)); lons=NaN(1,length(profiles));           % positions

    for i=1:length(profiles)
        % Copy stacked profiles
        nitrate_interp_all(:,i)=profiles(i).nitrate;
        chlorophyll_interp_all(:,i)=profiles(i).chlorophyll;

        % Integrate nitrate over depth (valid samples only)
        np=profiles(i).nitrate;
        valid=~isnan(np(depth_idx));
        if sum(valid)>=2
            d=common_depth(depth_idx(valid)); n=np(depth_idx(valid));
            integrated_nitrate(i)=trapz(d,n);                             % units: µmol/kg · m (pre-density)
        end

        % Integrate chlorophyll over depth
        cp=profiles(i).chlorophyll;
        valid_chla=~isnan(cp(depth_idx));
        if sum(valid_chla)>=2
            d=common_depth(depth_idx(valid_chla)); c=cp(depth_idx(valid_chla));
            integrated_chla(i)=trapz(d,c);                                % mg/m^2 if c in mg/m^3
        end

        % Time and position arrays
        dates_all(i)=profiles(i).date;
        lats(i)=profiles(i).latitude; lons(i)=profiles(i).longitude;
    end

    % Count how many profiles for this float had valid chlorophyll integrals
    chlorophyll_valid_counts(float_id)=sum(~isnan(integrated_chla));

    % Store full evolution for diagnostics (heatmaps & overlays)
    float_evolution(end+1)=struct('float_id',float_id,'dates',dates_all,'nitrate_profiles',nitrate_interp_all,'chlorophyll_profiles',chlorophyll_interp_all,'integrated_nitrate',integrated_nitrate,'integrated_chla',integrated_chla);

    % ============================================================
    % Per-float time series plot: integrated nitrate & chlorophyll
    % (only for floats listed in 'selected_floats')
    % ============================================================
    if any(strcmp(float_id,selected_floats))
        figure; hold on; yyaxis left;
        h1 = plot(dates_all, integrated_nitrate, '-k', 'LineWidth', 1.2); % integrated nitrate
        nitrate_years = year(dates_all);
        minHandles = []; maxHandles = [];

        % Mark seasonal windows (max/min) per year on the time series
        for y = unique(nitrate_years)
            idx_curr = nitrate_years==y; idx_prev = nitrate_years==y-1;
            max_idx = (month(dates_all)==12 & idx_prev) | (ismember(month(dates_all),1:3) & idx_curr);
            min_idx = ismember(month(dates_all),8:11) & idx_curr;

            % Plot max marker if present
            if any(max_idx)
                [~, i_max] = max(integrated_nitrate.*double(max_idx));
                maxHandles(end+1) = plot(dates_all(i_max), integrated_nitrate(i_max), 'r^', 'MarkerSize', 7, 'LineWidth', 1.2);
            end
            % Plot min marker if present
            if any(min_idx)
                [~, i_min] = min(integrated_nitrate.*double(min_idx) + (~min_idx)*1e10);
                minHandles(end+1) = plot(dates_all(i_min), integrated_nitrate(i_min), 'bo', 'MarkerSize', 7, 'LineWidth', 1.2);
            end
        end
        ylabel('Integrated Nitrate to 150 m (µ mol m^{-2})'); % label kept as in source

        % Chlorophyll on right axis (styled as green line with no markers)
        yyaxis right;
        hChl = plot(dates_all, integrated_chla, 'Color', [0 1 0], 'LineStyle', '-', 'LineWidth', 1.5, 'Marker', 'none', 'DisplayName', 'Integrated Chlorophyll');
        ylabel('Integrated Chlorophyll (mg m^{-2})'); xlabel('Date');
        title(sprintf('Integrated Nitrate and Chlorophyll - Float %s', float_id));

        % Construct legend robustly depending on which markers exist
        if ~isempty(minHandles) && ~isempty(maxHandles)
            legend([h1, minHandles(1), maxHandles(1), hChl], {'Integrated Nitrate','Min Nitrate','Max Nitrate','Integrated Chlorophyll'}, 'Location','best');
        else
            legend([h1, hChl], {'Integrated Nitrate','Integrated Chlorophyll'}, 'Location','best');
        end
        grid on;
    end

    % ================================================
    % Per-year ANCP computation for this float (robust)
    % ================================================
    nitrate_years=year(dates_all);
    unique_years=unique(nitrate_years);
    for y=unique_years
        idx_curr=nitrate_years==y; idx_prev=nitrate_years==y-1;

        % Seasonal windows
        max_idx=(month(dates_all)==12 & idx_prev)|(ismember(month(dates_all),1:3)&idx_curr);
        min_idx=ismember(month(dates_all),8:11)&idx_curr;

        % Require at least one candidate in both windows
        if sum(max_idx)<1 || sum(min_idx)<1, continue; end

        % NaN-robust selection of extreme indices in the windows
        maskMax=-inf(size(integrated_nitrate)); maskMax(max_idx)=integrated_nitrate(max_idx);
        [~,i_max]=max(maskMax,[], 'omitnan');
        if isempty(i_max)||isnan(integrated_nitrate(i_max)), continue; end

        maskMin=inf(size(integrated_nitrate)); maskMin(min_idx)=integrated_nitrate(min_idx);
        [~,i_min]=min(maskMin,[], 'omitnan');
        if isempty(i_min)||isnan(integrated_nitrate(i_min)), continue; end

        % Extract gridded profiles and compute profile-wise drawdown
        np_max=nitrate_interp_all(:,i_max);
        np_min=nitrate_interp_all(:,i_min);
        delta_nitrate=np_max-np_min;

        % Integrate drawdown over depth (valid samples only)
        valid=isfinite(delta_nitrate(depth_idx));
        if sum(valid)<2, continue; end
        d=common_depth(depth_idx(valid)); dn=delta_nitrate(depth_idx(valid));

        % ==== UNIT CONVERSION TO gC m^-2 ====
        % dn is in µmol/kg; multiply by seawater density (?) to get µmol/m^3,
        % convert to mol/m^3 (÷1e6), integrate over meters (mol/m^2),
        % then ×12 g/mol (molar mass of C) and ×R (C:N) to get gC/m^2.
        rho=1025; molar_mass_C=12;
        dn_molm3=dn.*rho/1e6;                 % µmol/kg × kg/m^3 / 1e6 = mol/m^3
        ANCP=trapz(d,dn_molm3)*molar_mass_C*R; % gC/m^2 (annual drawdown proxy)

        % Use the minimum-nitrate profile position for map placement (as requested)
        all_lat=[all_lat,lats(i_min)];
        all_lon=[all_lon,lons(i_min)];
        all_ANCP=[all_ANCP,ANCP]; all_years=[all_years,y]; float_ids=[float_ids,{float_id}];

        % Save the track (index span) for possible trajectory display
        ANCP_profile_lats{end+1}=lats(min(i_min,i_max):max(i_min,i_max));
        ANCP_profile_lons{end+1}=lons(min(i_min,i_max):max(i_min,i_max));

        % Store integrated chlorophyll at the minimum-nitrate time (for scatter)
        cp_min=chlorophyll_interp_all(:,i_min);
        valid_chla_min=~isnan(cp_min(depth_idx));
        if sum(valid_chla_min)>=2
            d_ch=common_depth(depth_idx(valid_chla_min));
            c_ch=cp_min(depth_idx(valid_chla_min));
            int_chla_min=trapz(d_ch,c_ch);
        else
            int_chla_min=NaN;
        end
        all_integrated_chla=[all_integrated_chla,int_chla_min];
    end
end

% =================
% === Global plots
% =================

% Scatter: ANCP vs integrated chlorophyll (linear)
figure; scatter(all_integrated_chla,all_ANCP,60,'filled');
xlabel('Integrated Chlorophyll (mg m^{-2})'); ylabel('ANCP (gC m^{-2} yr^{-1})');
title('ANCP vs Integrated Chlorophyll (Linear Scale)'); grid on;

% Scatter: ANCP vs integrated chlorophyll (log–log)
figure; loglog(all_integrated_chla,all_ANCP,'o','MarkerSize',6,'MarkerFaceColor','auto');
xlabel('Integrated Chlorophyll (mg m^{-2})'); ylabel('ANCP (gC m^{-2} yr^{-1})');
title('ANCP vs Integrated Chlorophyll (Log-Log Scale)'); grid on;

% Map: ANCP colored by value at the chosen minimum-nitrate locations
figure; geoscatter(all_lat,all_lon,60,all_ANCP,'filled'); colorbar;
title('Map of Annual Net Community Production (ANCP)'); colormap(jet); geobasemap colorterrain;

% Histogram: ANCP distribution with annotated literature ranges (illustrative bands)
figure;
histogram(all_ANCP, 'BinWidth', 10, 'FaceColor', [0.2 0.6 0.9], 'EdgeColor', 'k', 'LineWidth', 1.2);
xlabel('ANCP (gC m^{-2} yr^{-1})', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Histogram of ANCP Values', 'FontSize', 14, 'FontWeight', 'bold');
grid on; box on; hold on;
yl=ylim;
% Example bands (labels kept as in the source script)
plot([22.3 95.2],[yl(2)-2 yl(2)-2],'r-','LineWidth',3);
text(mean([22.3 95.2]),yl(2)-2.5,'Norwegian Sea (22.3–95.2)', ...
    'HorizontalAlignment','center','Color','r','FontSize',10,'FontWeight','bold');
plot([5 30],[yl(2)-5 yl(2)-5],'g-','LineWidth',3);
text(mean([5 30]),yl(2)-5.5,'Siberian Shelf (5–30)', ...
    'HorizontalAlignment','center','Color','g','FontSize',10,'FontWeight','bold');
hold off;

% Bar chart: number of valid chlorophyll integrals per float
keys_chl=keys(chlorophyll_valid_counts);
vals_chl=values(chlorophyll_valid_counts);
counts=cell2mat(vals_chl);
figure; bar(categorical(keys_chl),counts);
xlabel('Float ID'); ylabel('Valid Chlorophyll Profiles');
title('Number of Valid Chlorophyll Profiles per Float'); grid on;

% ===============================
% Profile evolution (heatmap) QA
% ===============================

% Nitrate evolution (up to 10 floats)
for i=1:min(10,length(float_evolution))
    evo=float_evolution(i);
    if isempty(evo.dates)||isempty(evo.nitrate_profiles), continue; end
    figure; imagesc(datenum(evo.dates),common_depth,evo.nitrate_profiles);
    set(gca,'YDir','reverse'); datetick('x','yyyy','keeplimits');
    xlabel('Year'); ylabel('Depth (m)');
    title(['Nitrate Profile Evolution - Float ' evo.float_id]);
    c=colorbar; c.Label.String='Nitrate (\mu mol/kg)'; % colorbar label (escaped µ in TeX)
end

% Chlorophyll evolution (up to 10 floats)
for i=1:min(10,length(float_evolution))
    evo=float_evolution(i);
    if isempty(evo.dates)||isempty(evo.chlorophyll_profiles), continue; end
    figure; imagesc(datenum(evo.dates),common_depth,evo.chlorophyll_profiles);
    set(gca,'YDir','reverse'); datetick('x','yyyy','keeplimits');
    xlabel('Year'); ylabel('Depth (m)');
    title(['Chlorophyll Profile Evolution - Float ' evo.float_id]);
    c=colorbar; c.Label.String='Chlorophyll (mg/m^3)';
end

% Final status message
fprintf('All floats processed. Total ANCP values calculated: %d\n',numel(all_ANCP));
