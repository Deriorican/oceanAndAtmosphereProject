from utils import *


data_folder = "C:\\Users\\lovat\\Desktop\\MyFiles\\Unif\\oc&at\\oceanAndAtmosphereProject\\data"
file_name = 'sample_20221201_1200.nc'
map_name = "ne_110m_admin_0_countries.shp"
save_folder = "C:\\Users\\lovat\\Desktop\\MyFiles\\Unif\\oc&at\\oceanAndAtmosphereProject\\plots"

file_path = pjoin(data_folder, file_name)
map_file = pjoin(data_folder, map_name)
ds = nc.Dataset(file_path)
print(ds.variables.keys())
geopot, geopot_dep = getVariableAndDependencies(ds, "z")



#500mbar ----------------------------------------------------------------------------------------------------------
#geopotential
savename_500mbar = "500mbar_geopotential.html"
savefile_500mbar = pjoin(save_folder, savename_500mbar)
plotOnSphere(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 0, :, :], map_file=map_file, title="500 mbar Geopotential", legend="Geopotential [m^2 s^-2]", save_file=savefile_500mbar)

#gradient intensity
savename_500mbar_grad = "500mbar_geopotential_grad.html"
savefile_500mbar_grad = pjoin(save_folder, savename_500mbar_grad)
grad_lat, grad_lon = computeGradient(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 0, :, :])
grad_lat = grad_lat.T
grad_lon = grad_lon.T
grad_norm = np.sqrt(grad_lat ** 2 + grad_lon ** 2)
plotOnSphere(geopot_dep['latitude'], geopot_dep["longitude"], grad_norm, map_file=map_file, title="500 mbar Geopotential Gradient", legend="Geopotential Gradient [m s^-2]", save_file=savefile_500mbar_grad)

#gradient as arrows + geopotential
savename_500mbar_grad_arrow = "500mbar_geopotential_grad_arrow.html"
savefile_500mbar_grad_arrow = pjoin(save_folder, savename_500mbar_grad_arrow)
dxdydz = prepareArrowForPlot(geopot_dep['latitude'], geopot_dep["longitude"], grad_lat, grad_lon)
plotOnSphere(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 0, :, :], map_file=map_file, title="500 mbar Geopotential with gradient", legend="Geopotential [m^2 s^-2]", save_file=savefile_500mbar_grad_arrow, arrows=dxdydz, partial_arrow_data=2e-2)



#1000mbar ---------------------------------------------------------------------------------------------------------
#geopotential
savename_1000mbar = "1000mbar_geopotential.html"
savefile_1000mbar = pjoin(save_folder, savename_1000mbar)
plotOnSphere(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 1, :, :], map_file=map_file, title="1000 mbar Geopotential", legend="Geopotential [m^2 s^-2]", save_file=savefile_1000mbar)

#gradient intensity
savename_1000mbar_grad = "1000mbar_geopotential_grad.html"
savefile_1000mbar_grad = pjoin(save_folder, savename_1000mbar_grad)
grad_lat, grad_lon = computeGradient(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 1, :, :])
grad_lat = grad_lat.T
grad_lon = grad_lon.T
grad_norm = np.sqrt(grad_lat ** 2 + grad_lon ** 2)
plotOnSphere(geopot_dep['latitude'], geopot_dep["longitude"], grad_norm, map_file=map_file, title="1000 mbar Geopotential Gradient", legend="Geopotential Gradient [m s^-2]", save_file=savefile_1000mbar_grad)

#gradient as arrows + geopotential
savename_1000mbar_grad_arrow = "1000mbar_geopotential_grad_arrow.html"
savefile_1000mbar_grad_arrow = pjoin(save_folder, savename_1000mbar_grad_arrow)
dxdydz = prepareArrowForPlot(geopot_dep['latitude'], geopot_dep["longitude"], grad_lat, grad_lon)
plotOnSphere(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 1, :, :], map_file=map_file, title="1000 mbar Geopotential with gradient", legend="Geopotential [m^2 s^-2]", save_file=savefile_1000mbar_grad_arrow, arrows=dxdydz, partial_arrow_data=2e-2)
