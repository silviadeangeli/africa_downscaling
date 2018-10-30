from rasterizza_analizza import analyse
gar= "/Users/silvia/Documents/AFRICA_DATA/Equatorial_Guinea/GAR_original/africa_gnq_only_buiA.shp"
builtup = "/Users/silvia/Documents/AFRICA_DATA/Equatorial_Guinea/GQ_buiA_GHSL_GUF_LC_OSM_wp3_90m.tif"
results_tif= "/Users/silvia/Documents/AFRICA_DATA/Equatorial_Guinea/count.tif"
valfis_tif= "/Users/silvia/Documents/AFRICA_DATA/Equatorial_Guinea/count_valfis.tif"
results_txt= "/Users/silvia/Documents/AFRICA_DATA/Equatorial_Guinea/count.txt"

print(analyse(gar, builtup, results_tif, valfis_tif, results_txt))