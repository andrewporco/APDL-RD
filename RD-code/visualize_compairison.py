import pyvista as pv
import numpy as np

# Load the two .npz files
data1 = np.load('Part3_material.npz')
data2 = np.load('Part3_netfabb.npz')

# Function to create surface from vertices and displacements
def create_surface(verts, disp, additional_scale=1, name=""):
    scale_factor = 1 * additional_scale
    verts_scaled = verts * scale_factor
    disp_scaled = disp * scale_factor
    
    disp_magnitude = np.linalg.norm(disp_scaled, axis=1)
    
    # Normalize displacement magnitudes to range [0, 1]
    disp_magnitude_normalized =disp_magnitude
    
    print(f"\n{name} Displacement Magnitude Stats:")
    if name == 'Netfabb':
        disp_magnitude *= 1000
    print(f"Min: {np.min(disp_magnitude):.6f}")
    print(f"Max: {np.max(disp_magnitude):.6f}")
    print(f"Mean: {np.mean(disp_magnitude):.6f}")
    print(f"Std Dev: {np.std(disp_magnitude):.6f}")
    
    cloud = pv.PolyData(verts_scaled)
    cloud['Displacement Magnitude'] = disp_magnitude
    cloud['Normalized Displacement'] = disp_magnitude_normalized
    
    surface = cloud.delaunay_3d().extract_surface()
    surface_with_data = surface.sample(cloud)
    
    return surface_with_data

# Create surfaces for both datasets
surface1 = create_surface(data1['verts'], data1['disp'], additional_scale=1, name="APDL")
surface2 = create_surface(data2['verts'], data2['disp'], additional_scale=0.0025, name="Netfabb")

# Create a plotter with two subplots side by side
plotter = pv.Plotter(shape=(1, 2))

# Function to add mesh and center it
def add_centered_mesh(plotter, surface, subplot_index, title):
    plotter.subplot(0, subplot_index)
    
    plotter.add_mesh(surface, scalars='Normalized Displacement', cmap='jet', 
                     show_edges=False, smooth_shading=True)  # Set color limits to [0, 1] for normalized data
    
    plotter.add_scalar_bar(title='Normalized Displacement',
                           label_font_size=10,
                           width=0.25,
                           position_x=.75)
    
    # Calculate the range manually
    disp_mag = surface['Displacement Magnitude']
    actual_range = (np.min(disp_mag), np.max(disp_mag))
    
    plotter.add_text(f"{title}\nActual Range: {actual_range[0]:.2e} - {actual_range[1]:.2e} mm",
                     position='upper_left', font_size=10)
    
    # Center the camera on the mesh
    plotter.camera_position = 'iso'
    plotter.camera.focal_point = surface.center
    plotter.camera.position = surface.center + np.array([0, 0, surface.length])
    plotter.camera.up = (0, 1, 0)
    plotter.reset_camera()

# Add the first shape to the left subplot
add_centered_mesh(plotter, surface1, 0, "APDL")

# Add the second shape to the right subplot
add_centered_mesh(plotter, surface2, 1, "Netfabb")

# Show the plot
plotter.show()