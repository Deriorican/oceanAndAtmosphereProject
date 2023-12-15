from utils import *
import matplotlib.pyplot as plt


data_folder = "C:\\Users\\lovat\\Desktop\\MyFiles\\Unif\\oc&at\\oceanAndAtmosphereProject\\data"
file_name = 'sample_20221201_1200.nc'
map_name = "ne_110m_admin_0_countries.shp"
save_folder = "C:\\Users\\lovat\\Desktop\\MyFiles\\Unif\\oc&at\\oceanAndAtmosphereProject\\plots"
savename_500mbar = "500mbar_geostrophic.html"
savename_1000mbar = "1000mbar_geostrophic.html"

savename_500mbar_th = "500mbar_geostrophic_theoretical.html"
savename_1000mbar_th = "1000mbar_geostrophic_theoretical.html"

savename_500mbar_diff = "500mbar_geostrophic_theoreticalAndMeasureDiff.html"
savename_1000mbar_diff = "1000mbar_geostrophic_theoreticalAndMeasureDiff.html"

file_path = pjoin(data_folder, file_name)
map_file = pjoin(data_folder, map_name)
savefile_500mbar_th = pjoin(save_folder, savename_500mbar_th)
savefile_1000mbar_th = pjoin(save_folder, savename_1000mbar_th)
savefile_500mbar = pjoin(save_folder, savename_500mbar)
savefile_1000mbar = pjoin(save_folder, savename_1000mbar)
savefile_500mbar_diff = pjoin(save_folder, savename_500mbar_diff)
savefile_1000mbar_diff = pjoin(save_folder, savename_1000mbar_diff)

ds = nc.Dataset(file_path)
print(ds.variables.keys())
geopot, geopot_dep = getVariableAndDependencies(ds, "z")
percentage = 0.02


#COMPARISON
geowind_zonal, geowind_zonal_dep = getVariableAndDependencies(ds, "u")
geowind_meridional, geowind_meridional_dep = getVariableAndDependencies(ds, "v")

#comparison 500
grad_lat, grad_lon = computeGradient(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 0, :, :])
grad_lat = grad_lat.T; grad_lon = grad_lon.T
grad_norm = np.sqrt(grad_lat**2 + grad_lon**2)
grad = np.array([grad_lon, grad_lat])
rad_lat = geopot_dep['latitude'] / 90 * np.pi/2
rad_lon = geopot_dep["longitude"] / 180 * np.pi
lat_mesh, _ = np.meshgrid(rad_lat, rad_lon)
lat_mesh = lat_mesh.T
geostrophic_wind = computeGeostrophicWind(lat_mesh, grad)
mean_measured_wind = np.mean(np.sqrt(geowind_zonal[0, 0, :, :]**2 + geowind_meridional[0, 0, :, :]**2))
mean_theoretical_wind = np.mean(np.sqrt(geostrophic_wind[0]**2 + geostrophic_wind[1]**2))
delta_angle = 180 / np.pi * np.arccos((geostrophic_wind[0] * geowind_zonal[0, 0, :, :] + geostrophic_wind[1] * geowind_meridional[0, 0, :, :]) / (np.sqrt(geostrophic_wind[0]**2 + geostrophic_wind[1]**2) * np.sqrt(geowind_zonal[0, 0, :, :]**2 + geowind_meridional[0, 0, :, :]**2)))
delta_norm = np.abs((np.sqrt(geostrophic_wind[0]**2 + geostrophic_wind[1]**2) - np.sqrt(geowind_zonal[0, 0, :, :]**2 + geowind_meridional[0, 0, :, :]**2))) / mean_measured_wind

RMSE_angle = np.sqrt(np.mean(delta_angle**2))
RMSE_norm = np.sqrt(np.mean(delta_norm**2))
print(f"500 mbar :\n\tRMSE angle : {RMSE_angle}\n\tRMSE norm : {RMSE_norm}")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 15))
fig.suptitle('Geostrophic wind differences between theoretical and measured data (500mbar)')
flat_lat = lat_mesh.flatten()
flat_grad_norm = grad_norm.flatten()
flat_delta = delta_angle.flatten()
flat_delta_norm = delta_norm.flatten()
selected_indices = np.random.choice(len(flat_lat), size = int(percentage * len(flat_lat)), replace=False)
flat_lat = flat_lat[selected_indices] * 180 / np.pi
flat_grad_norm = flat_grad_norm[selected_indices]
flat_delta = flat_delta[selected_indices]
flat_delta_norm = flat_delta_norm[selected_indices]
                  
ax1.scatter(flat_lat[np.abs(flat_lat) >= 20], flat_delta[np.abs(flat_lat) >= 20], alpha = 0.2, label = "valid latitudes")
ax2.scatter(flat_grad_norm[np.abs(flat_lat) >= 20], flat_delta[np.abs(flat_lat) >= 20], alpha = 0.2, label = "valid latitudes")
ax1.scatter(flat_lat[np.abs(flat_lat) < 20], flat_delta[np.abs(flat_lat) < 20], alpha = 0.2, color="red", label = "unvalid latitudes")
ax2.scatter(flat_grad_norm[np.abs(flat_lat) < 20], flat_delta[np.abs(flat_lat) < 20], alpha = 0.2, color="red", label = "unvalid latitudes")

ax3.scatter(flat_lat[np.abs(flat_lat) >= 20], flat_delta_norm[np.abs(flat_lat) >= 20], alpha = 0.2, label = "valid latitudes")
ax4.scatter(flat_grad_norm[np.abs(flat_lat) >= 20], flat_delta_norm[np.abs(flat_lat) >= 20], alpha = 0.2, label = "valid latitudes")
ax3.scatter(flat_lat[np.abs(flat_lat) < 20], flat_delta_norm[np.abs(flat_lat) < 20], alpha = 0.2, color="red", label = "unvalid latitudes")
ax4.scatter(flat_grad_norm[np.abs(flat_lat) < 20], flat_delta_norm[np.abs(flat_lat) < 20], alpha = 0.2, color="red", label = "unvalid latitudes")

ax1.set_title("-Variations with latitude-")
ax2.set_title("-Variations with geopotential gradient-")
ax1.set_xlabel("latitude [rad]")
ax2.set_xlabel("geopotential gradient norm [m/rad]")
ax1.set_ylabel("Angular difference of the wind direction [°]")
ax2.set_ylabel("Angular difference of the wind direction [°]")

ax3.set_title("-Variations with latitude-")
ax4.set_title("-Variations with geopotential gradient-")
ax3.set_xlabel("latitude [rad]")
ax4.set_xlabel("geopotential gradient norm [m/rad]")
ax3.set_ylabel("Relative difference in the norm\nof the geostrophic wind [-]")
ax4.set_ylabel("Relative difference in the norm\nof the geostrophic wind [-]")

ax1.grid(linestyle="--")
ax2.grid(linestyle="--")
ax3.grid(linestyle="--")
ax4.grid(linestyle="--")
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

fig.tight_layout(pad=5.0)
plt.show()

#comp 1000
grad_lat, grad_lon = computeGradient(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 1, :, :])
grad_lat = grad_lat.T; grad_lon = grad_lon.T
grad_norm = np.sqrt(grad_lat**2 + grad_lon**2)
grad = np.array([grad_lon, grad_lat])
rad_lat = geopot_dep['latitude'] / 90 * np.pi/2
rad_lon = geopot_dep["longitude"] / 180 * np.pi
lat_mesh, _ = np.meshgrid(rad_lat, rad_lon)
lat_mesh = lat_mesh.T
geostrophic_wind = computeGeostrophicWind(lat_mesh, grad)
mean_measured_wind = np.mean(np.sqrt(geowind_zonal[0, 1, :, :]**2 + geowind_meridional[0, 1, :, :]**2))
mean_theoretical_wind = np.mean(np.sqrt(geostrophic_wind[0]**2 + geostrophic_wind[1]**2))
print(mean_measured_wind, mean_theoretical_wind)
delta_angle = 180 / np.pi * np.arccos((geostrophic_wind[0] * geowind_zonal[0, 1, :, :] + geostrophic_wind[1] * geowind_meridional[0, 1, :, :]) / (np.sqrt(geostrophic_wind[0]**2 + geostrophic_wind[1]**2) * np.sqrt(geowind_zonal[0, 1, :, :]**2 + geowind_meridional[0, 1, :, :]**2)))
delta_norm = np.abs((np.sqrt(geostrophic_wind[0]**2 + geostrophic_wind[1]**2) - np.sqrt(geowind_zonal[0, 1, :, :]**2 + geowind_meridional[0, 1, :, :]**2))) / mean_measured_wind

RMSE_angle = np.sqrt(np.mean(delta_angle**2))
RMSE_norm = np.sqrt(np.mean(delta_norm**2))
print(f"1000 mbar :\n\tRMSE angle : {RMSE_angle}\n\tRMSE norm : {RMSE_norm}")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 15))
fig.suptitle('Geostrophic wind differences between theoretical and measured data (1000mbar)')
flat_lat = lat_mesh.flatten()
flat_grad_norm = grad_norm.flatten()
flat_delta = delta_angle.flatten()
flat_delta_norm = delta_norm.flatten()
selected_indices = np.random.choice(len(flat_lat), size = int(percentage * len(flat_lat)), replace=False)
flat_lat = flat_lat[selected_indices] * 180 / np.pi
flat_grad_norm = flat_grad_norm[selected_indices]
flat_delta = flat_delta[selected_indices]
flat_delta_norm = flat_delta_norm[selected_indices]
                  
ax1.scatter(flat_lat[np.abs(flat_lat) >= 20], flat_delta[np.abs(flat_lat) >= 20], alpha = 0.2, label = "valid latitudes")
ax2.scatter(flat_grad_norm[np.abs(flat_lat) >= 20], flat_delta[np.abs(flat_lat) >= 20], alpha = 0.2, label = "valid latitudes")
ax1.scatter(flat_lat[np.abs(flat_lat) < 20], flat_delta[np.abs(flat_lat) < 20], alpha = 0.2, color="red", label = "unvalid latitudes")
ax2.scatter(flat_grad_norm[np.abs(flat_lat) < 20], flat_delta[np.abs(flat_lat) < 20], alpha = 0.2, color="red", label = "unvalid latitudes")

ax3.scatter(flat_lat[np.abs(flat_lat) >= 20], flat_delta_norm[np.abs(flat_lat) >= 20], alpha = 0.2, label = "valid latitudes")
ax4.scatter(flat_grad_norm[np.abs(flat_lat) >= 20], flat_delta_norm[np.abs(flat_lat) >= 20], alpha = 0.2, label = "valid latitudes")
ax3.scatter(flat_lat[np.abs(flat_lat) < 20], flat_delta_norm[np.abs(flat_lat) < 20], alpha = 0.2, color="red", label = "unvalid latitudes")
ax4.scatter(flat_grad_norm[np.abs(flat_lat) < 20], flat_delta_norm[np.abs(flat_lat) < 20], alpha = 0.2, color="red", label = "unvalid latitudes")

ax1.set_title("-Variations with latitude-")
ax2.set_title("-Variations with geopotential gradient-")
ax1.set_xlabel("latitude [rad]")
ax2.set_xlabel("geopotential gradient norm [m/rad]")
ax1.set_ylabel("Angular difference of the wind direction [°]")
ax2.set_ylabel("Angular difference of the wind direction [°]")

ax3.set_title("-Variations with latitude-")
ax4.set_title("-Variations with geopotential gradient-")
ax3.set_xlabel("latitude [rad]")
ax4.set_xlabel("geopotential gradient norm [m/rad]")
ax3.set_ylabel("Relative difference in the norm\nof the geostrophic wind [-]")
ax4.set_ylabel("Relative difference in the norm\nof the geostrophic wind [-]")

ax1.grid(linestyle="--")
ax2.grid(linestyle="--")
ax3.grid(linestyle="--")
ax4.grid(linestyle="--")
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

fig.tight_layout(pad=5.0)
plt.show()


#theoretical 500mbar
grad_lat, grad_lon = computeGradient(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 0, :, :])
grad_lat = grad_lat.T; grad_lon = grad_lon.T
grad = np.array([grad_lon, grad_lat])
rad_lat = geopot_dep['latitude'] / 90 * np.pi/2
rad_lon = geopot_dep["longitude"] / 180 * np.pi
lat_mesh, _ = np.meshgrid(rad_lat, rad_lon)
lat_mesh = lat_mesh.T
geostrophic_wind = computeGeostrophicWind(lat_mesh, grad)
geostrophic_wind_norm = np.sqrt(geostrophic_wind[0]**2 + geostrophic_wind[1]**2)
dxdydz = dxdydz = prepareArrowForPlot(geopot_dep['latitude'], geopot_dep["longitude"], geostrophic_wind[1]*1e-3, geostrophic_wind[0]*1e-3)
plotOnSphere(geopot_dep['latitude'], geopot_dep["longitude"], geostrophic_wind_norm, map_file=map_file, title="500 mbar Geostrophic wind", legend="Geostrophic wind m/s", save_file=savefile_500mbar_th, arrows=dxdydz, partial_arrow_data=2e-2, scaleref=1.0)

#theoretical 1000mbar
grad_lat, grad_lon = computeGradient(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 1, :, :])
grad_lat = grad_lat.T; grad_lon = grad_lon.T
grad = np.array([grad_lon, grad_lat])
rad_lat = geopot_dep['latitude'] / 90 * np.pi/2
rad_lon = geopot_dep["longitude"] / 180 * np.pi
lat_mesh, _ = np.meshgrid(rad_lat, rad_lon)
lat_mesh = lat_mesh.T
geostrophic_wind = computeGeostrophicWind(lat_mesh, grad)
geostrophic_wind_norm = np.sqrt(geostrophic_wind[0]**2 + geostrophic_wind[1]**2)
dxdydz = dxdydz = prepareArrowForPlot(geopot_dep['latitude'], geopot_dep["longitude"], geostrophic_wind[1]*1e-3, geostrophic_wind[0]*1e-3)
plotOnSphere(geopot_dep['latitude'], geopot_dep["longitude"], geostrophic_wind_norm, map_file=map_file, title="1000 mbar Geostrophic wind", legend="Geostrophic wind m/s", save_file=savefile_1000mbar_th, arrows=dxdydz, partial_arrow_data=2e-2, scaleref=.5)

#measure
geowind_zonal, geowind_zonal_dep = getVariableAndDependencies(ds, "u")
geowind_meridional, geowind_meridional_dep = getVariableAndDependencies(ds, "v")

#500mbar
geostrophic_wind_norm = np.sqrt(geowind_zonal[0, 0, :, :]**2 + geowind_meridional[0, 0, :, :]**2)
dxdydz = dxdydz = prepareArrowForPlot(geowind_zonal_dep['latitude'], geowind_zonal_dep["longitude"], geowind_meridional[0, 0, :, :]*1e-3, geowind_zonal[0, 0, :, :]*1e-3)
plotOnSphere(geowind_zonal_dep['latitude'], geowind_zonal_dep["longitude"], geostrophic_wind_norm, map_file=map_file, title="500 mbar Geostrophic wind", legend="Geostrophic wind m/s", save_file=savefile_500mbar, arrows=dxdydz, partial_arrow_data=2e-2, scaleref=1.0)

#1000mbar
geostrophic_wind_norm = np.sqrt(geowind_zonal[0, 1, :, :]**2 + geowind_meridional[0, 1, :, :]**2)
dxdydz = dxdydz = prepareArrowForPlot(geowind_zonal_dep['latitude'], geowind_zonal_dep["longitude"], geowind_meridional[0, 1, :, :]*1e-3, geowind_zonal[0, 1, :, :]*1e-3)
plotOnSphere(geowind_zonal_dep['latitude'], geowind_zonal_dep["longitude"], geostrophic_wind_norm, map_file=map_file, title="1000 mbar Geostrophic wind", legend="Geostrophic wind m/s", save_file=savefile_1000mbar, arrows=dxdydz, partial_arrow_data=2e-2, scaleref=1.5)


#diff 500mbar
grad_lat, grad_lon = computeGradient(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 0, :, :])
grad_lat = grad_lat.T; grad_lon = grad_lon.T
grad = np.array([grad_lon, grad_lat])
rad_lat = geopot_dep['latitude'] / 90 * np.pi/2
rad_lon = geopot_dep["longitude"] / 180 * np.pi
lat_mesh, _ = np.meshgrid(rad_lat, rad_lon)
lat_mesh = lat_mesh.T
geostrophic_wind = computeGeostrophicWind(lat_mesh, grad)
diff_norm = np.sqrt((geostrophic_wind[0] - geowind_zonal[0, 0, :, :])**2 + (geostrophic_wind[1]-geowind_meridional[0, 0, :, :])**2)
plotOnSphere(geowind_zonal_dep['latitude'], geowind_zonal_dep["longitude"], diff_norm, map_file=map_file, title="1000 mbar Geostrophic wind differences\nBetween measures and theory", legend="Norm of the difference [m/s]", save_file=savefile_500mbar_diff)

#diff 1000mbar
grad_lat, grad_lon = computeGradient(geopot_dep['latitude'], geopot_dep["longitude"], geopot[0, 1, :, :])
grad_lat = grad_lat.T; grad_lon = grad_lon.T
grad = np.array([grad_lon, grad_lat])
rad_lat = geopot_dep['latitude'] / 90 * np.pi/2
rad_lon = geopot_dep["longitude"] / 180 * np.pi
lat_mesh, _ = np.meshgrid(rad_lat, rad_lon)
lat_mesh = lat_mesh.T
geostrophic_wind = computeGeostrophicWind(lat_mesh, grad)
diff_norm = np.sqrt((geostrophic_wind[0] - geowind_zonal[0, 1, :, :])**2 + (geostrophic_wind[1]-geowind_meridional[0, 1, :, :])**2)
plotOnSphere(geowind_zonal_dep['latitude'], geowind_zonal_dep["longitude"], diff_norm, map_file=map_file, title="1000 mbar Geostrophic wind differences\nBetween measures and theory", legend="Norm of the difference [m/s]", save_file=savefile_1000mbar_diff)