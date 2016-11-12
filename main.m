%author: NasaFireProject

%%%%%%%%%%%%%%NOTES%%%%%%%%%%%%%%%%%%%%%
%How it works:
%for a given region of the planet, look for:
%latlim: latitude limits
%monlim: longitude limits
%years: period of investigation (must be included in 1997-2015)
%geozone: description of the zone under study. 'indonesia' for example.
%wet_months: list of the months of the wet season

%create a function analysis_{geozone you are investigating} calling
%the function geozone_analysis with the proper arguments (cf below)

%add this function to the main() function

%run the script and enjoy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function main()
%     analysis_indonesia();
    analysis_australia();
    %...
end

%complete analysis on different regions of the globe
function analysis_indonesia()
    geozone_analysis([-10 10], [94 144], [1997, 2015], 'indonesia', {'01','02','03','04','10','11','12'});
end

function analysis_australia()
    geozone_analysis([-44 -9], [112 154], [1997, 2015], 'australia', {'01','02','03','04','10','11','12'});
end

%complete analysis of the fire activity on a given region
%latlim: latitude limits
%lonlim: longitude limits
%years: years limit
%geozone: description of the zone under study. 'indonesia' for example.
%wet_months: list of months of the wet season
function geozone_analysis(latlim, lonlim, years, geozone, wet_months)
    import map.*
    cols = [721+4*lonlim(1), 720+4*lonlim(2)];
    rows = [361+4*latlim(1), 360+4*latlim(2)];
    [bigmat, row_names] = create_bigmat(years, rows, cols);  %#ok<*ASGLU>
    bigmat_norm = normalize_bigmat(bigmat, row_names);
   
    %charts
    monthly_total_norm = monthly_total_ba(bigmat_norm);
    plot_monthly_total_ba(monthly_total_norm, years, geozone);
    
    monthly_total = monthly_total_ba(bigmat);
    monthly_mean = monthly_mean_ba(monthly_total, row_names, years);
    plot_annual_cycle(monthly_mean,geozone);
    
    [compo1, compo2, compo3, pca_timeseries] = perform_pca(bigmat, rows, cols);
    plot_pca_timeseries(pca_timeseries, years);
    
    %maps
    [wet_season_ba, dry_season_ba] = get_seasonal_ba(row_names, bigmat, wet_months, rows, cols);
    
        %plot 12 monthly maps
    plot_monthly_maps(row_names, bigmat, latlim, lonlim, rows, cols, geozone)
    
        %plot 2 seasonal maps
    label='Burned area (m^2)';
    plot_map(wet_season_ba, latlim, lonlim, strcat(geozone, '_wetseason'), label);
    plot_map(dry_season_ba, latlim, lonlim, strcat(geozone, '_dryseason'), label);
       
        %plot 3 principal composant maps
    label='';
    plot_map(compo1, latlim, lonlim, strcat(geozone, '_PC1'),label);
    plot_map(compo2, latlim, lonlim, strcat(geozone, '_PC2'),label);
    plot_map(compo3, latlim, lonlim, strcat(geozone, '_PC3'),label);
end

%Aggregation of HDF5 datasets
function [bigmat, row_names] = create_bigmat(years, rows, cols)
    nb_months = (years(2) - years(1) + 1)*12;
    nb_locs = (cols(2)-cols(1)+1)*(rows(2)-rows(1)+1);
    bigmat = zeros(nb_months,nb_locs);
    cpt = 1;
    row_names = {};
    for year = years(1):years(2)
        file_path = strcat('data/GFED4.1s_',num2str(year),'.hdf5');
        dataset_name_surface = char(strcat('/ancill/grid_cell_area'));
        data_surface = h5read(file_path, dataset_name_surface);
        for month = {'01','02','03','04','05','06','07','08','09','10','11','12'}
            dataset_name_frac = char(strcat('/burned_area/',month,'/burned_fraction'));
            data_frac = h5read(file_path, dataset_name_frac);
            data = data_frac.*data_surface;
            data = data';
            data = flipud(data);
            indo_matrix = data(rows(1):rows(2),cols(1):cols(2));
            row_indo_matrix = reshape(indo_matrix, 1, nb_locs);
            bigmat(cpt,:) = row_indo_matrix;
            row_names = [row_names, strcat(month, '_', num2str(year))];  %#ok<*AGROW>
            cpt = cpt + 1;
        end
    end
end

%normalization of the array: retrieving monthly means for each location to
%each element of the array
function bigmat_norm = normalize_bigmat(bigmat, row_names)

    bigmat_norm = bigmat;
    dims = size(bigmat);
    nb_locs = dims(2);
    nb_months = dims(1);
    avgs = zeros(12, nb_locs);
    %Building the monthly mean matrix (12xnb_locs)
    for k = 1:nb_months
        name = row_names(k);
        nam = char(name);
        nam_spl = strsplit(nam, '_');
        month_spl = char(nam_spl(1));
        month = str2double(month_spl);
        avgs(month,:) = avgs(month,:) + bigmat(k,:);
    end
    avgs = avgs / 19;
    %Retrieving monthly means for each location
    for k = 1:nb_months
        name = row_names(k);
        nam = char(name);
        nam_spl = strsplit(nam, '_');
        month_spl = char(nam_spl(1));
        month = str2double(month_spl);
        bigmat_norm(k,:) = bigmat_norm(k,:) - avgs(month,:);
    end
end

%Seasons: Create two matrices: one for the dry season, one for wet season
function [wet_season_ba, dry_season_ba] = get_seasonal_ba(row_names, mat, wet_months, rows, cols)
    dims = size(mat);
    nb_locs = dims(2);
    nb_months = dims(1);
    nb_rows = rows(2) - rows(1) + 1;
    nb_cols = cols(2) - cols(1) + 1;    
    
    wet_bigmat = zeros(1, nb_locs);
    dry_bigmat = zeros(1, nb_locs);
    for k = 1:nb_months
        name = row_names(k);
        nam = char(name);
        arr = strsplit(nam, '_');

        if ismember(arr(1),wet_months) == 1
            wet_bigmat(1,:) = wet_bigmat(1,:) + mat(k,:);
        else
            dry_bigmat(1,:) = dry_bigmat(1,:) + mat(k,:);
        end  
    end

    wet_season_ba = reshape(wet_bigmat, nb_rows, nb_cols);
    dry_season_ba = reshape(dry_bigmat, nb_rows, nb_cols);
    
    %csvwrite('csv_files/dry_season_ba.csv',dry_season_ba);
    %csvwrite('csv_files/wet_season_ba.csv',wet_season_ba);
end

%Create matrix nb_monthsx1 with burned areas per month through the years
function y = monthly_total_ba(mat)
    y = sum(mat,2);
end

%Create matrix 12x1 with total burned area per month
function y = monthly_mean_ba(monthly_total, row_names, years)
    nb_years = years(2) - years(1) + 1;
    y = zeros(12,1);
    for k = 1:(nb_years*12)
        name = row_names(k);
        nam = char(name);
        nam_spl = strsplit(nam, '_');
        month_spl = char(nam_spl(1));
        month = str2double(month_spl);
        y(month) = y(month) + monthly_total(k);
    end
    y = y/nb_years;
end

%Plots
function plot_annual_cycle(monthly_mean_ba, geozone) 
    fig=figure;
    % Define values for x, y1, and y2
    x1 = 1:12;
    y1 = monthly_mean_ba;

    % Plot y1 vs. x (blue, solid) and y2 vs. x (red, dashed)
    
    plot(x1, y1)

    % Turn on the grid
    grid on

    % Set the axis limits
    axis([1 12 -inf inf])
    
    %Set x axis names
    set(gca, 'XTick', 1:12, 'XTickLabel', {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'})
    
    % Add title and axis labels
    titlo=['Annual cycle - average monthly burned area - ' geozone];
    title(titlo)
    xlabel('Month')
    ylabel('Total burned area (m^2)')
    saveas(fig,['Output/' titlo '.png'])
end

function plot_monthly_total_ba(monthly_total_ba, years, geozone)
    fig=figure;
    nb_months = size(monthly_total_ba);
    nb_months = nb_months(1);
    % Define values for x, y1, and y2
    
    x1 = 1:nb_months;
    x1 = x1/12;
    x1 = x1+years(1);
    y1 = monthly_total_ba;

    elnino = csvread('data/elnino.txt');
    elnino_ts = reshape(elnino', 1, nb_months);
    [hAx,hLine1,hLine2] = plotyy(x1, y1, x1, elnino_ts);
    
    % Turn on the grid
    grid on

    % Set the axis limits
    axis([1997 2016 -inf inf])
    
    
    % Add title and axis labels
    titlo=['Deviation of monthly burned area' geozone];
    title(titlo);
    xlabel('Month')
    ylabel(hAx(1),'Total burned area (m^2)')
    ylabel(hAx(2),'Oceanic Niño Index (ONI)')
    saveas(fig,['Output/' titlo '.png'])
end

%Perform PCA
function [compo1, compo2, compo3, pca_timeseries] = perform_pca(bigmat, rows, cols) %#ok<*DEFNU>

    coeff = pca(bigmat);
    compo123 = coeff(:,1:3);
    pca_timeseries = bigmat * compo123;
    disp(pca_timeseries);
    nb_rows = rows(2) - rows(1) + 1;
    nb_cols = cols(2) - cols(1) + 1;
   
    compo1 = reshape(coeff(:,1), nb_rows, nb_cols);
    compo2 = reshape(coeff(:,2), nb_rows, nb_cols);
    compo3 = reshape(coeff(:,3), nb_rows, nb_cols);

end

%plot the 3 first PC and their timeseries
function plot_pca_timeseries(pca_timeseries, years)
    fig=figure;
    nb_months = 12*(years(2)-years(1)+1);
 % Define values for x, y1, and y2
    x = 1:nb_months;
    x = x/12;
    x = x+years(1);
    y1 = pca_timeseries(:,1);
    y2 = pca_timeseries(:,2);
    y3 = pca_timeseries(:,3);

    % Plot y1 vs. x (blue, solid) and y2 vs. x (red, dashed)
%     figure
    plot(x, y1, x, y2, x, y3)
    

    % Turn on the grid
    grid on

    % Set the axis limits
    axis([years(1) years(2)+1 -inf inf])


    % Add title and axis labels
    titlo = char('Time series of Principal Component 1, 2 & 3');
    title(titlo)
    xlabel('Month')
    ylabel('Burned area (m^2)')
    saveas(fig,['Output/' titlo '.png'])
end

%(unused)plot the elnino trend on the period 1997-2015
function plot_elnino()
    fig=figure;
    elnino = csvread('data/elnino.txt');
    elnino_ts = reshape(elnino', 1, 228);    

    x = 1:228;
    x = x/12;
    x = x+1997;
    y0 = zeros(1,228);
    
    plot(x, elnino_ts, x, y0)
    grid on
    axis([1997 2016 -2.5 4])
    
    titlo = char('ElNino Index monthly value 1997 - 2015');
    title(titlo)
    xlabel('Month')
    ylabel('Oceanic Niño Index (ONI)')
    saveas(fig,['Output/' titlo '.png'])
end

%creates the map for a given matrix and the limts of lat and lon
function fig=plot_map(matrix, latlim, lonlim, titlo, label)
    fig=figure('Color','white');
    dims = size(matrix);
    disp(dims)
%     figure('Color','white')
    title(titlo);
    [lat, lon]= meshgrat(latlim, lonlim, [dims(1), dims(2)]);

    worldmap(latlim, lonlim);
    load coastlines
    surfacem(lat, lon, matrix);
    c = colorbar('southoutside');
    c.Label.String = label;
    plotm(coastlat, coastlon);
    saveas(fig,['Output/' titlo '.png'])
end

function plot_monthly_maps(row_names, mat, latlim, lonlim, rows, cols, geozone)
    nb_rows = rows(2) - rows(1) + 1;
    nb_cols = cols(2) - cols(1) + 1;
    dims = size(mat);
    nb_locs = dims(2);
    nb_months = dims(1);
    for k = 1:12
        monthly_mat = zeros(1,nb_locs);
        for j = 1:nb_months
            name = row_names(j);
            nam = char(name);
            arr = strsplit(nam, '_');
            if str2double(arr(1)) == k
                monthly_mat(1,:) = monthly_mat(1,:) + mat(j,:);
            end
        end
        monthly_mat = reshape(monthly_mat, nb_rows, nb_cols);
        titre = strcat('Month #', num2str(k), ' -  ', geozone);
        label='Burned area (m^2)';
        plot_map(monthly_mat, latlim, lonlim, titre, label);
    end
    
end

function annualMatrix = findAnnualAvg(final_matrix, years)
    dims = size(final_matrix);
    nb_years = years(2)-years(1)+1;
    annualMatrix=sum(final_matrix, 2)/nb_years;
    annualMatrix=reshape(annualMatrix, dims(1), dims(2));
end

function monthlyMatrices=calcMonthlyAvg(final_matrix, years)
    monthlyMatrices=cell(12,1);
    nb_years = years(2)-years(1)+1;

    for i=1:12
        tempMonthMat=sum(final_matrix(:, i:12:end), 2)/nb_years;
        monthlyMatrices{i}=tempMonthMat;
    end
end

