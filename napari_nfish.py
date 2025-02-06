import napari
from magicgui import magicgui, widgets
import nd2
from skimage import filters, measure
from magicgui.widgets import FileEdit
import os
import numpy as np
from scipy.ndimage import maximum_filter
from cellpose import models
import torch
import pandas as pd

# Spot detection function
def analyze_nfish_spots(image_data, threshImageSpot=0.002, previous_segmentation=None):
    if previous_segmentation is not None:
        segmented_DAPI = previous_segmentation
    else:
        segmented_DAPI = cellposeSegNucG(image_data[0])
    
    spots_detected1 = detect_spots(image_data[1], threshImageSpot)
    spots_df = extract_spot_info(spots_detected1, image_data[1], segmented_DAPI)
    df_nuc_inten = calculate_intensity_per_label(image_data[1], segmented_DAPI)

    # Count cells and spots
    cell_count = np.unique(segmented_DAPI).size - 1  # Exclude background
    spot_count = len(spots_df)
    
    # Calculate spots per cell
    spots_per_cell = spot_count / cell_count if cell_count > 0 else 0

    # Print results
    print(f"Cell Count: {cell_count}")
    print(f"Spot Count: {spot_count}")
    print(f"Spots per Cell: {spots_per_cell:.2f}")
    print(f"threshold used for this spot detection: {threshImageSpot}")

    return spots_df, df_nuc_inten, segmented_DAPI, spots_detected1, cell_count, spot_count

# Detect spots in the FISH image
def detect_spots(image_to_analyze, threshold):
    sigma = 2.0
    log_image = filters.laplace(filters.gaussian(image_to_analyze, sigma=sigma))
    max_filtered_image = maximum_filter(log_image, size=3)
    spots_detected = (max_filtered_image > threshold)
    return spots_detected

# Segmentation of nuclei using Cellpose
def cellposeSegNucG(img, min_size=15):
    use_gpu = torch.cuda.is_available()
    model = models.Cellpose(gpu=use_gpu, model_type="nuclei")
    res, _, _, _ = model.eval(img, channels=[0, 0], diameter=None, min_size=min_size)
    return res

# Extract information about each FISH spot
def extract_spot_info(binary_spots_image, nuclei_image, segmented_nuclei_image):
    spots_info = []
    labeled_spots, _ = measure.label(binary_spots_image, return_num=True, connectivity=1, background=0)
    properties = measure.regionprops(labeled_spots, intensity_image=nuclei_image)
    
    for prop in properties:
        y, x = prop.centroid
        intensity = prop.mean_intensity
        cell_label = segmented_nuclei_image[int(y), int(x)]
        spots_info.append({'Y': y, 'X': x, 'Intensity': intensity, 'CellLabel': cell_label})
    
    spots_df = pd.DataFrame(spots_info)
    return spots_df

# Mock function for intensity calculation
def calculate_intensity_per_label(image, segmented_image):
    props = measure.regionprops_table(
        measure.label(segmented_image),
        intensity_image=image,
        properties=['label', 'mean_intensity']
    )
    df = pd.DataFrame(props)
    return df

# Create the main viewer
viewer = napari.Viewer()
previous_segmentation = None  # Variable to store the previous segmentation

# GUI to load nd2 files with a file dialog box
@magicgui(call_button="Load Image", layout='vertical')
def load_button_widget():
    global previous_segmentation
    previous_segmentation = None  # Reset previous segmentation when a new image is loaded
    file_edit = FileEdit(mode="r", filter="*.nd2")
    file_edit.show()
    file_edit.changed.connect(lambda: load_file(file_edit.value))

# Function to load and display the ND2 file
def load_file(filepath):
    if filepath and os.path.isfile(filepath):
        # Remove existing layers
        while len(viewer.layers) > 0:
            viewer.layers.pop()
        
        im_read = nd2.imread(filepath)  # Load the ND2 file
        viewer.add_image(im_read, name='image')  # Add image to the viewer

# Slider for threshold
threshold_slider = widgets.FloatSlider(name="Threshold", min=0.0001, max=0.005, step=0.0001, value=0.001)
# GUI to detect spots with adjustable threshold

from qtpy.QtWidgets import QLabel
label_widget = QLabel()
label_widget.setText("")
label_widget.setStyleSheet("color: white; font-size: 16px;")
viewer.window.add_dock_widget(label_widget, name="Spot and Cell Counts", area='bottom')

# GUI to detect spots with adjustable threshold
@magicgui(call_button="Detect Spots", layout='vertical')
def detect_button_widget():
    global previous_segmentation

    if 'image' not in viewer.layers:
        print("No image loaded. Please load an image first.")
        return

    # Remove existing "Detected Spots" layer if it exists
    if 'Detected Spots' in viewer.layers:
        viewer.layers.remove('Detected Spots')

    image_layer = viewer.layers['image']
    image_data = image_layer.data

    # Get threshold from slider
    threshold = threshold_slider.value

    # Use previous segmentation if available
    spots_df, df_nuc_inten, imSeg, spots_detected2, total_cell, total_spot = analyze_nfish_spots(
        image_data, threshold, previous_segmentation)

    # Store the current segmentation for future use
    if previous_segmentation is None:
        previous_segmentation = imSeg

    # Only add the "Detected Spots" layer
    viewer.add_labels(spots_detected2, name='Detected Spots')

    # Display total spot and cell count in the Napari viewer
    label_widget.setText(f"Total Cells: {total_cell}\nTotal Spots: {total_spot}\nSpots/Cell: {total_spot/total_cell if total_cell > 0 else 0:.2f}")

# Add widgets to Napari
viewer.window.add_dock_widget(load_button_widget, name="Load Image", area='right')
viewer.window.add_dock_widget(threshold_slider, name="Threshold", area='right')
viewer.window.add_dock_widget(detect_button_widget, name="Detect Spots", area='right')

# Start Napari viewer
napari.run()
