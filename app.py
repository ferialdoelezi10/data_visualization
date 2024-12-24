from flask import Flask, render_template
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
from netCDF4 import Dataset
import geopandas as gpd
import os

app = Flask(__name__)

@app.route("/")
def index():
    # Part 1: Line Chart for CO₂ Concentrations Over Time
    data = pd.read_excel("co2_data.xlsx", header=0, names=['Year', 'CO2'])
    data = data[(data['Year'] >= 1959) & (data['Year'] <= 2023)]  # Filter data if needed

    # Create Plotly figure for line chart
    line_fig = go.Figure()
    line_fig.add_trace(go.Scatter(
        x=data['Year'],
        y=data['CO2'],
        mode='lines+markers',
        line=dict(color='blue'),
        marker=dict(size=8, color='blue', symbol='circle'),
        name='CO2 Concentration'
    ))

    # Add animation frames
    frames = [
        go.Frame(
            data=[
                go.Scatter(
                    x=data['Year'][:i],
                    y=data['CO2'][:i],
                    mode='lines+markers',
                    line=dict(color='blue'),
                    marker=dict(size=8, color='blue', symbol='circle')
                )
            ],
            name=str(data['Year'].iloc[i - 1])
        )
        for i in range(1, len(data) + 1)
    ]
    line_fig.frames = frames

    # Configure play/pause buttons and slider
    line_fig.update_layout(
        updatemenus=[
            dict(
                type='buttons',
                showactive=True,
                buttons=[
                    dict(
                        label="Play",
                        method="animate",
                        args=[None, dict(frame=dict(duration=700, redraw=True), fromcurrent=True)]
                    ),
                    dict(
                        label="Pause",
                        method="animate",
                        args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")]
                    )
                ]
            )
        ],
        sliders=[
            dict(
                active=0,
                steps=[
                    dict(
                        method="animate",
                        args=[
                            [str(data['Year'].iloc[k])],
                            dict(mode="immediate", frame=dict(duration=500, redraw=True), transition=dict(duration=0))
                        ],
                        label=str(data['Year'].iloc[k])
                    )
                    for k in range(len(data))
                ],
                x=0, y=0,
                currentvalue=dict(font=dict(size=12), prefix="Year: ", visible=True),
                len=1.0
            )
        ]
    )

    # Layout for line chart
    line_fig.update_layout(
        title="Average carbon dioxide (CO₂) levels in the atmosphere worldwide from 1959 to 2023",
        xaxis=dict(title="Year", tickmode='linear', dtick=5, range=[1959, 2023]),
        yaxis=dict(title="CO2 Concentration (ppm)", range=[310, 420]),
        hovermode="x unified"
    )

    # Convert line chart to HTML
    line_graph_html = pio.to_html(line_fig, full_html=False)

    # Part 2: 3D Globe Visualization
    # Paths to NetCDF and GeoJSON files
    nc_file_path = "/Users/aldoelezi/Desktop/3D_concentration_2022/3D_concentration_2022.nc"
    geojson_file_path = "/Users/aldoelezi/Desktop/adm.geojson"

    # Load NetCDF data
    data_nc = Dataset(nc_file_path, mode='r')
    lats = data_nc.variables['latitude'][:]
    lons = data_nc.variables['longitude'][:]
    co2_levels = data_nc.variables['co2_conc'][0, 0, :, :]  # Adjust based on your file

    # Fix longitude wrapping
    lons[lons >= 180] -= 360

    # Create latitude and longitude grid
    lon, lat = np.meshgrid(lons, lats)

    # Close the longitude gap to avoid the line in the middle
    lon = np.hstack([lon, lon[:, 0:1] + 360])
    lat = np.hstack([lat, lat[:, 0:1]])
    co2_levels = np.hstack([co2_levels, co2_levels[:, 0:1]])

    # Map latitude and longitude to 3D spherical coordinates
    R = 1
    x = R * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
    y = R * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
    z = R * np.sin(np.radians(lat))

    # Load and process GeoJSON
    gdf = gpd.read_file(geojson_file_path)
    if gdf.crs != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")

    gdf['geometry'] = gdf.simplify(tolerance=0.001)
    gdf['geometry'] = gdf.boundary

    # Generate border traces
    border_traces = []
    for _, row in gdf.iterrows():
        geometry = row.geometry

        if geometry.geom_type == 'LineString':
            coords = list(geometry.coords)
            lon_border, lat_border = zip(*coords)
            lon_border = np.array(lon_border)
            lon_border[lon_border >= 180] -= 360
            x_border = R * np.cos(np.radians(lat_border)) * np.cos(np.radians(lon_border))
            y_border = R * np.cos(np.radians(lat_border)) * np.sin(np.radians(lon_border))
            z_border = R * np.sin(np.radians(lat_border))
            border_traces.append(go.Scatter3d(
                x=x_border,
                y=y_border,
                z=z_border,
                mode='lines',
                line=dict(color='black', width=1)
            ))

        elif geometry.geom_type == 'MultiLineString':
            for line in geometry.geoms:
                coords = list(line.coords)
                lon_border, lat_border = zip(*coords)
                lon_border = np.array(lon_border)
                lon_border[lon_border >= 180] -= 360
                x_border = R * np.cos(np.radians(lat_border)) * np.cos(np.radians(lon_border))
                y_border = R * np.cos(np.radians(lat_border)) * np.sin(np.radians(lon_border))
                z_border = R * np.sin(np.radians(lat_border))
                border_traces.append(go.Scatter3d(
                    x=x_border,
                    y=y_border,
                    z=z_border,
                    mode='lines',
                    line=dict(color='black', width=1)
                ))

    # Add CO₂ surface with hover text
    globe_fig = go.Figure(go.Surface(
        x=x, y=y, z=z,
        surfacecolor=co2_levels,
        colorscale=[[0, "blue"], [0.3, "yellow"], [0.4, "orange"], [1, "red"]],
        cmin=405.31384, cmax=440.4423,
        showscale=True,
        colorbar=dict(title="CO₂ ppm")
    ))

    # Add border traces to globe
    for trace in border_traces:
        globe_fig.add_trace(trace)

    # Layout for 3D globe
    globe_fig.update_layout(
        title="3D Earth with CO₂ Concentrations and Borders",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        ),
        showlegend=False  # Disable the legend
    )

    # Convert 3D globe to HTML
    globe_graph_html = pio.to_html(globe_fig, full_html=False)

    # Render the HTML template with both graphs
    return render_template(
        'index.html',
        line_graph_html=line_graph_html,
        globe_graph_html=globe_graph_html
    )

if __name__ == "__main__":
    app.run(debug=True)
