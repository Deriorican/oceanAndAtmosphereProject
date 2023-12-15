import numpy as np
import pandas as pd
import geopandas as gpd
import netCDF4 as nc
from os.path import join as pjoin
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def computeCurl(lat, lon, data_lat, data_lon):
    R = 6371000.0
    lat_ex = np.zeros((len(lat) + 2, len(lon) + 2))
    lon_ex = np.zeros((len(lat) + 2, len(lon) + 2))
    data_lat_ex = np.zeros((len(lat) + 2, len(lon) + 2))
    data_lon_ex = np.zeros((len(lat) + 2, len(lon) + 2))
    lon_ex[1:-1, 1:-1], lat_ex[1:-1, 1:-1] = np.meshgrid(lon * np.pi / 180, lat * np.pi / 180)
    data_lat_ex[1:-1, 1:-1] = data_lat
    data_lon_ex[1:-1, 1:-1] = data_lon
    lat_ex[0, :] = np.pi - np.roll(lat_ex[2, :], len(lon) // 2); lat_ex[-1, :] = -np.pi - np.roll(lat_ex[-3, :], len(lon) // 2)
    lon_ex[:, 0] = -(2 * np.pi - lon_ex[:, -2]); lon_ex[:, -1] = 2 * np.pi + lon_ex[:, 1]
    data_lat_ex[0, :] = np.roll(data_lat_ex[2, :], len(lon) // 2); data_lat_ex[-1, :] = np.roll(data_lat_ex[-3, :], len(lon) // 2)
    data_lat_ex[:, 0] = data_lat_ex[:, -2]; data_lat_ex[:, -1] = data_lat_ex[:, 1]
    data_lon_ex[0, :] = np.roll(data_lon_ex[2, :], len(lon) // 2); data_lon_ex[-1, :] = np.roll(data_lon_ex[-3, :], len(lon) // 2)
    data_lon_ex[:, 0] = data_lon_ex[:, -2]; data_lon_ex[:, -1] = data_lon_ex[:, 1]
    lat_comp = (data_lon_ex[0:-2, 1:-1] * np.cos(lat_ex[0:-2, 1:-1]) - data_lon_ex[2:, 1:-1] * np.cos(lat_ex[2:, 1:-1])) / (lat_ex[0:-2, 1:-1] - lat_ex[2:, 1:-1])
    lon_comp = np.zeros(data_lat.shape)
    lon_comp[1:-1, :] = (data_lat_ex[2:-2, 2:] - data_lat_ex[2:-2, :-2]) / (lon_ex[2:-2, 2:] - lon_ex[2:-2, :-2])
    res = (1 / (R * np.cos(lat_ex[1:-1, 1:-1]))) * (lat_comp - lon_comp)
    return res




def computeGeostrophicWind(lat, geopotential_grad, lat_threshold=20):
    f = 2 * 7.2921159e-5 * np.sin(lat)
    threshold_rad = lat_threshold * np.pi / 180
    f = np.where(np.abs(lat) < threshold_rad, 1.0, f)
    u0 = np.where(np.abs(lat) < threshold_rad, 0.0, - (1 / f) * geopotential_grad[1])
    u1 = np.where(np.abs(lat) < threshold_rad, 0.0, (1 / f) * geopotential_grad[0])
    return  np.array([u0, u1])




def computeGradient(lat, lon, data):
    R = 6371000.0
    lat_ex = np.zeros((len(lat) + 2, len(lon) + 2))
    lon_ex = np.zeros((len(lat) + 2, len(lon) + 2))
    data_ex = np.zeros((len(lat) + 2, len(lon) + 2))
    lon_ex[1:-1, 1:-1], lat_ex[1:-1, 1:-1] = np.meshgrid(lon * np.pi / 180, lat * np.pi / 180)
    data_ex[1:-1, 1:-1] = data
    lat_ex[0, :] = np.pi - np.roll(lat_ex[2, :], len(lon) // 2); lat_ex[-1, :] = -np.pi - np.roll(lat_ex[-3, :], len(lon) // 2)
    lon_ex[:, 0] = -(2 * np.pi - lon_ex[:, -2]); lon_ex[:, -1] = 2 * np.pi + lon_ex[:, 1]
    data_ex[0, :] = np.roll(data_ex[2, :], len(lon) // 2); data_ex[-1, :] = np.roll(data_ex[-3, :], len(lon) // 2)
    data_ex[:, 0] = data_ex[:, -2]; data_ex[:, -1] = data_ex[:, 1]
    dlat_data = (data_ex[0:-2, 1:-1] - data_ex[2:, 1:-1]) / (R * (lat_ex[0:-2, 1:-1] - lat_ex[2:, 1:-1]))
    dlon_data = np.zeros(data.shape)
    dlon_data[1:-1, :] = (data_ex[2:-2, 2:] - data_ex[2:-2, :-2]) / (R * (lon_ex[2:-2, 2:] - lon_ex[2:-2, :-2]) * np.cos(lat_ex[2:-2, 1:-1]))
    return dlat_data.T, dlon_data.T


def prepareArrowForPlot(lat, lon, data_lat, data_lon):
    rad_lat = lat / 90 * np.pi/2
    rad_lon = lon / 180 * np.pi
    u, v = np.meshgrid(rad_lat, rad_lon)
    u = u.T; v = v.T
    x = np.cos(u) * np.cos(v)
    y = np.cos(u) * np.sin(v)
    z = np.sin(u)
    u_ = u + data_lat
    v_ = v + data_lon
    x_ = np.cos(u_) * np.cos(v_)
    y_ = np.cos(u_) * np.sin(v_)
    z_ = np.sin(u_)
    dx = x_ - x
    dy = y_ - y
    dz = z_ - z
    return dx, dy, dz



def getVariableAndDependencies(ds, variable_name):
    variable_values = ds.variables[variable_name][:]
    dimensions = ds.variables[variable_name].dimensions
    dependencies = {}
    for dimension in dimensions:
        dependencies[dimension] =  ds.variables[dimension][:]
    return variable_values, dependencies



def plot_polygon(poly):
    
    xy_coords = poly.exterior.coords.xy
    lon = np.array(xy_coords[0])
    lat = np.array(xy_coords[1])
    
    lon = lon * np.pi/180
    lat = lat * np.pi/180
    
    R = 1 + 1e-3
    x = R* np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)
    
    return x, y, z



def plotOnSphere(lat, lon, data, title=" ", legend=" ", map_file=None, save_file=None, arrows=None, html=True, partial_arrow_data=1.0, colorscale=None, scaleref=1.0):            
    rad_lat = lat / 90 * np.pi/2
    rad_lon = lon / 180 * np.pi
    u, v = np.meshgrid(rad_lat, rad_lon)
    u = u.T; v = v.T
    x = np.cos(u) * np.cos(v)
    y = np.cos(u) * np.sin(v)
    z = np.sin(u)
    fig = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    subplot_titles=[title],
                    )
    if colorscale is None:
        fig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=data, colorbar_x=0), 1, 1)
    else :
        fig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=data, colorbar_x=0, cmin=colorscale[0], cmax=colorscale[1]), 1, 1)

    if arrows is not None:
        dx = arrows[0][25:-25, :]
        dy = arrows[1][25:-25, :]
        dz = arrows[2][25:-25, :]
        x = x[25:-25, :]
        y = y[25:-25, :]
        z = z[25:-25, :]
        selected_indices = np.random.choice(len(x.flatten()), size = int(partial_arrow_data * len(x.flatten())), replace=False)
        print(f"{len(selected_indices)} vectors displayed.")
        fig.add_trace(go.Cone(x=x.flatten()[selected_indices] * 1.01, y=y.flatten()[selected_indices] * 1.01, z=z.flatten()[selected_indices] * 1.01,\
                              u=dx.flatten()[selected_indices] * 1.01, v=dy.flatten()[selected_indices] * 1.01, w=dz.flatten()[selected_indices] * 1.01,\
                              sizemode="scaled", sizeref=scaleref, anchor="center", colorscale=[[0, 'black'], [1, 'black']]), 1, 1)

    if map_file is not None:
        gdf = gpd.read_file(map_file)
        for i in gdf.index :
            polys = gdf.loc[i].geometry         # Polygons or MultiPolygons 
            if polys.geom_type == 'Polygon':
                x, y, z = plot_polygon(polys)
                fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color=f'rgb(0, 0,0)'), showlegend=False), 1, 1 ) 
            elif polys.geom_type == 'MultiPolygon':
                for poly in polys.geoms:
                    x, y, z = plot_polygon(poly)
                    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color=f'rgb(0, 0,0)'), showlegend=False), 1, 1 ) 
        

    fig.update_layout(title_text=legend)
    if save_file is not None:
        if html:
            fig.write_html(save_file)
        else : 
            fig.write_image(save_file)
    fig.show()
