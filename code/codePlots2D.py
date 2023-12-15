import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def compute_geopotential_gradient(lon, lat, geopotential):
    gradient_lon, gradient_lat = np.gradient(geopotential, axis=(1, 0))
    return gradient_lon, gradient_lat

def calculate_geostrophic_wind(geopotential, dx, dy, f):
    d_geopotential_dx, d_geopotential_dy = np.gradient(geopotential, dx, dy)
    v_geo = (-1/f * d_geopotential_dy)/rayon_Terre
    u_geo = (1/f * d_geopotential_dx)/rayon_Terre
    return u_geo, v_geo
rayon_Terre = 6371000.0
# Load the NetCDF file
file_path = 'C:/Users/Thomas/Desktop/atmo/p2/sample_20221201_1200.nc'
data = xr.open_dataset(file_path)

# Extract data at 500 hPa level and a specific time (adjust time accordingly)
geopotential_500hPa = data['z'].sel(level=500, time='2022-12-01T12:00:00')

# Subsample data for better performance
subsample_factor = 25
lon_sub = geopotential_500hPa.longitude[::subsample_factor]
lat_sub = geopotential_500hPa.latitude[::subsample_factor]
geopotential_sub = geopotential_500hPa.values[::subsample_factor, ::subsample_factor]

# Reverse the latitude array
#lat_sub = np.flipud(lat_sub)
#lon_sub = np.flipud(lon_sub)
# Create a grid of longitude and latitude directly from the xarray dataset
lon, lat = np.meshgrid(lon_sub, lat_sub)

# Compute the gradient of the geopotential
grad_geopotential_x, grad_geopotential_y = compute_geopotential_gradient(lon, lat, geopotential_sub)

# Calculate geostrophic wind
dx = np.radians(0.25)
dy = np.radians(0.25)
f = 2 * 7.2921e-5 * np.sin(np.radians(lat))  # Coriolis parameter

u_geo, v_geo = calculate_geostrophic_wind(geopotential_sub, dx, dy, f)

# Plotting
fig1 = plt.figure(figsize=(10, 7))
ax1 = fig1.add_subplot(111, projection=ccrs.PlateCarree())
contour_geopotential = geopotential_500hPa.plot.contourf(ax=ax1, transform=ccrs.PlateCarree(), cmap='viridis', levels=20)
ax1.coastlines()
ax1.gridlines()
ax1.set_title('Geopotential at 500 hPa - 2022-12-01 12:00 UTC')

fig2 = plt.figure(figsize=(10, 7))
ax2 = fig2.add_subplot(111, projection=ccrs.PlateCarree())
quiver_geopotential = ax2.quiver(lon, lat, grad_geopotential_x, grad_geopotential_y, scale=40000, color='red', transform=ccrs.PlateCarree())
ax2.coastlines()
ax2.gridlines()
ax2.set_title('Geopotential Gradient (500 hPa) Quiver Plot - 2022-12-01 12:00 UTC')
# Calculate the magnitude of geostrophic wind
magnitude_geo_wind = np.sqrt(u_geo**2 + v_geo**2)

# Plotting
fig3 = plt.figure(figsize=(10, 7))
ax3 = fig3.add_subplot(111, projection=ccrs.PlateCarree())

# Filter data to exclude values between latitudes -20 and 20 degrees
mask = (lat < -10) | (lat > 10)
u_geo_filtered = np.where(mask, u_geo, np.nan)
v_geo_filtered = np.where(mask, v_geo, np.nan)
magnitude_geo_wind_filtered = np.where(mask, magnitude_geo_wind, np.nan)

# Plot quiver plot with colormap representing magnitude
quiver_geostrophic_wind = ax3.quiver(lon, lat, u_geo_filtered, v_geo_filtered, scale=8000, color='blue', transform=ccrs.PlateCarree())
contour_magnitude = ax3.contourf(lon, lat, magnitude_geo_wind_filtered, cmap='viridis', levels=20, alpha=0.5, transform=ccrs.PlateCarree())

ax3.coastlines()
ax3.gridlines()
ax3.set_title('Geostrophic Wind Magnitude (500 hPa) Quiver Plot - 2022-12-01 12:00 UTC')

# Add colorbar
cbar = plt.colorbar(contour_magnitude, ax=ax3, orientation='vertical', fraction=0.046, pad=0.04)
cbar.set_label('Magnitude of Geostrophic Wind (m/s)')

plt.show()



