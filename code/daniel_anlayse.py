from utils import *
import matplotlib.pyplot as plt


data_folder = "C:\\Users\\lovat\\Desktop\\MyFiles\\Unif\\oc&at\\oceanAndAtmosphereProject\\data"
file_name = 'sample_20230905_1200.nc'
map_name = "ne_110m_admin_0_countries.shp"
save_folder = "C:\\Users\\lovat\\Desktop\\MyFiles\\Unif\\oc&at\\oceanAndAtmosphereProject\\plots"
savename_daniel_500mbar_curl = "daniel_vorticity_500mbar.html"
savename_daniel_1000mbar_curl = "daniel_vorticity_1000mbar.html"

file_path = pjoin(data_folder, file_name)
map_file = pjoin(data_folder, map_name)
savefile_daniel_500mbar = pjoin(save_folder, savename_daniel_500mbar_curl)
savefile_daniel_1000mbar = pjoin(save_folder, savename_daniel_1000mbar_curl)

ds = nc.Dataset(file_path)
print(ds.variables.keys())

geowind_zonal, geowind_zonal_dep = getVariableAndDependencies(ds, "u")
geowind_meridional, geowind_meridional_dep = getVariableAndDependencies(ds, "v")

vorticity_500 = computeCurl(geowind_zonal_dep['latitude'], geowind_zonal_dep['longitude'], geowind_zonal[0, 0, :, :], geowind_meridional[0, 0, :, :])
vorticity_1000 = computeCurl(geowind_zonal_dep['latitude'], geowind_zonal_dep['longitude'], geowind_zonal[0, 1, :, :], geowind_meridional[0, 1, :, :])

lats, lons = np.meshgrid(geowind_zonal_dep['longitude'] * np.pi / 180, geowind_zonal_dep['latitude'] * np.pi / 180)


plotOnSphere(geowind_zonal_dep['latitude'], geowind_zonal_dep['longitude'], np.abs(vorticity_500), map_file=map_file, title="500 mbar vorticity", legend="vorticity [1/s]", save_file=savefile_daniel_500mbar, colorscale=[0, 2e-4])
plotOnSphere(geowind_zonal_dep['latitude'], geowind_zonal_dep['longitude'], np.abs(vorticity_1000), map_file=map_file, title="1000 mbar vorticity", legend="vorticity [1/s]", save_file=savefile_daniel_1000mbar, colorscale=[0, 1e-4])

