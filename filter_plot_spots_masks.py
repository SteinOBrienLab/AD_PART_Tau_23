import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def find_matching_mask(spot_filename, mask_files):
    """
    Find corresponding segmentation mask file for a given spot file.
    Args:
        spot_filename (str): Name of the spots file
        mask_files (list): List of mask file paths
    Returns:
        str: Path to matching mask file, or None if no match found
    """
    # Split filename into parts based on underscores
    parts = spot_filename.split('_')
    # Find index of marker identifier (AC)
    ac_index = parts.index('AC')  
    # Get image number that comes before AC
    image_number = parts[ac_index - 1]
    # Find part starting with 'orig' to get the original image identifier
    orig_part = next(part for part in parts if part.startswith('orig'))
    suffix = orig_part.replace('orig', '')
    
    # Construct expected mask filename pattern using image number and suffix
    mask_pattern = f"MAX_{image_number}_AC_NT{suffix}_seg.npy"
    print(f"Image number: {image_number}")
    print(f"Looking for mask matching pattern: {mask_pattern}")
    
    # Search for matching mask file
    matching_masks = [m for m in mask_files if os.path.basename(m) == mask_pattern]
    if matching_masks:
        return matching_masks[0]
    return None

def filter_spots_by_mask(csv_path, mask_path, output_path, spot_filename):
    """
    Filter spots based on cell segmentation mask and calculate per-cell statistics.
    Args:
        csv_path (str): Path to spots CSV file
        mask_path (str): Path to segmentation mask file
        output_path (str): Path to save filtered results
        spot_filename (str): Original spots filename for metadata
    Returns:
        tuple: (filtered spots DataFrame, ROI statistics list)
    """
    # Load spots data and mask
    spots = pd.read_csv(csv_path)
    mask = np.load(mask_path, allow_pickle=True).item()
    mask = mask['masks']
    
    # Get dimensions of mask image
    img_height, img_width = mask.shape
    
    # Get coordinate ranges of spots
    x_min, x_max = spots['x'].min(), spots['x'].max()
    y_min, y_max = spots['y'].min(), spots['y'].max()
    
    # Convert spot coordinates to pixel coordinates in mask space
    x_pixels = ((spots['x'] - x_min) / (x_max - x_min) * img_width).round().astype(int)
    y_pixels = ((spots['y'] - y_min) / (y_max - y_min) * img_height).round().astype(int)
    
    # Filter out spots outside image boundaries
    valid_coords = (x_pixels >= 0) & (x_pixels < mask.shape[1]) & \
                  (y_pixels >= 0) & (y_pixels < mask.shape[0])
    spots_filtered = spots[valid_coords].copy()
    x_pixels = x_pixels[valid_coords]
    y_pixels = y_pixels[valid_coords]
    
    # Get cell IDs for each spot location
    mask_values = mask[y_pixels, x_pixels]
    
    # Keep only spots inside cells (mask value > 0)
    spots_filtered = spots_filtered[mask_values > 0]
    x_pixels = x_pixels[mask_values > 0]
    y_pixels = y_pixels[mask_values > 0]
    
    # Add cell ID and pixel coordinates to filtered spots
    spots_filtered['cell_id'] = mask[y_pixels, x_pixels]
    spots_filtered['x_pixels'] = x_pixels
    spots_filtered['y_pixels'] = y_pixels
    
    # Get image metadata from filename
    parts = spot_filename.split('_')
    ac_index = parts.index('AC')
    image_number = parts[ac_index - 1]
    
    # Calculate median spot intensity for normalization
    median_spot_intensity = spots_filtered['intensity'].median()
    
    # Calculate statistics for each ROI/cell
    roi_stats = []
    unique_rois = np.unique(mask)
    unique_rois = unique_rois[unique_rois > 0]  # Exclude background (0)
    
    # Process each ROI
    for roi in unique_rois:
        # Calculate area of ROI
        roi_mask = mask == roi
        roi_area = np.sum(roi_mask)
        
        # Get spots in this ROI
        roi_spots = spots_filtered[spots_filtered['cell_id'] == roi]
        spots_in_roi = len(roi_spots)
        
        # Calculate intensity metrics, normalized count is generated by dividing the total ROI/neuron
        # intensity by the median spot intensity for that image
        total_intensity = roi_spots['intensity'].sum()
        normalized_count = total_intensity / median_spot_intensity if median_spot_intensity > 0 else 0
        
        # Get original image identifier
        orig_part = next(part for part in parts if part.startswith('orig'))
        
        # Store ROI statistics
        roi_stats.append({
            'Image_Number': image_number,
            'Orig_Number': f"{image_number}_AC_{orig_part}",
            'ROI_Key': f"{image_number}_AC_{orig_part}_ROI{roi}",
            'ROI_ID': roi,
            'Area_Pixels': roi_area,
            'APLP2': spots_in_roi,  # Gene name/marker count
            'Total_Intensity': total_intensity,
            'Median_Intensity': median_spot_intensity,
            'Normalized_Count': normalized_count,
            'Spots_per_1000px': (spots_in_roi / roi_area) * 1000
        })
    
    # Save filtered spots and statistics
    spots_filtered.to_csv(output_path, index=False)
    roi_stats_path = output_path.replace('filtered.csv', 'roi_stats.csv')
    pd.DataFrame(roi_stats).to_csv(roi_stats_path, index=False)
    
    return spots_filtered, roi_stats

def plot_spots_and_masks(spots, mask, filtered_spots, output_path, title="", point_size=3, show_plot=False):
    """
    Create visualization of spots overlaid on cell masks.
    Args:
        spots (DataFrame): Original spots data
        mask (ndarray): Segmentation mask
        filtered_spots (DataFrame): Filtered spots data
        output_path (str): Path to save plot
        title (str): Plot title
        point_size (int): Size of spot points
        show_plot (bool): Whether to display plot
    """
    # Get mask dimensions
    img_height, img_width = mask.shape
    
    # Convert spot coordinates to pixel space
    x_min, x_max = spots['x'].min(), spots['x'].max()
    y_min, y_max = spots['y'].min(), spots['y'].max()
    x_pixels = ((spots['x'] - x_min) / (x_max - x_min) * img_width).round().astype(int)
    y_pixels = ((spots['y'] - y_min) / (y_max - y_min) * img_height).round().astype(int)
    
    # Create plot
    plt.figure(figsize=(15, 12))
    
    # Plot cell masks with different colors
    plt.imshow(mask, cmap='tab20')
    
    # Plot all spots in red
    plt.scatter(x_pixels, y_pixels, c='red', s=point_size, alpha=0.3, label='All spots')
    
    # Plot filtered spots in yellow
    plt.scatter(filtered_spots['x_pixels'], filtered_spots['y_pixels'], 
               c='yellow', s=point_size, alpha=0.7, label='Filtered spots')
    
    # Add title and statistics
    plt.title(f'{title}\nTotal spots: {len(spots)}, Filtered spots: {len(filtered_spots)}\nNumber of cells: {len(np.unique(mask))-1}')
    plt.legend()
    plt.axis('off')
    
    # Save and optionally display plot
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    if show_plot:
        plt.show()
    plt.close()

def process_and_plot_spots(spots_dir, masks_dir, output_dir):
    """
    Process all spot files, match with masks, and create visualizations.
    Args:
        spots_dir (str): Directory containing spot files
        masks_dir (str): Directory containing mask files
        output_dir (str): Directory to save outputs
    """
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all spot and mask files
    spot_files = glob.glob(os.path.join(spots_dir, "RadialSymmetry_results_*.csv"))
    mask_files = glob.glob(os.path.join(masks_dir, "*_seg.npy"))
    
    print(f"Found {len(spot_files)} spot files and {len(mask_files)} mask files")
    
    # Initialize data collection
    summary_data = []
    combined_roi_stats = []
    
    # Process each spot file
    for spot_file in spot_files:
        spot_filename = os.path.basename(spot_file)
        print(f"\nProcessing spot file: {spot_filename}")
        
        # Find matching mask file
        mask_file = find_matching_mask(spot_filename, mask_files)
        if mask_file is None:
            print(f"No matching mask found for {spot_filename}")
            continue
            
        print(f"Found matching mask: {os.path.basename(mask_file)}")
        
        # Setup output paths
        base_name = spot_filename.split('.tif.csv')[0]
        filtered_csv = os.path.join(output_dir, f"{base_name}_filtered.csv")
        plot_path = os.path.join(output_dir, f"{base_name}_overlay.png")
        
        try:
            # Filter spots and generate statistics
            filtered_spots, roi_stats = filter_spots_by_mask(spot_file, mask_file, filtered_csv, spot_filename)
            combined_roi_stats.extend(roi_stats)
            
            # Load data for plotting
            spots = pd.read_csv(spot_file)
            mask = np.load(mask_file, allow_pickle=True).item()['masks']
            
            # Create visualization
            plot_spots_and_masks(spots, mask, filtered_spots, plot_path,
                               title=base_name, point_size=3, show_plot=False)
            
            # Extract image number for summary
            parts = spot_filename.split('_')
            ac_index = parts.index('AC')
            image_number = parts[ac_index - 1]
            
            # Calculate per-cell statistics
            spots_per_cell = filtered_spots.groupby('cell_id').size()
            
            # Add summary statistics
            summary_data.append({
                'Image': base_name,
                'Image_Number': image_number,
                'Total_Spots': len(spots),
                'Filtered_Spots': len(filtered_spots),
                'Cells_With_Spots': len(filtered_spots['cell_id'].unique()),
                'Max_Spots_Per_Cell': spots_per_cell.max() if not spots_per_cell.empty else 0,
                'Mean_Spots_Per_Cell': spots_per_cell.mean() if not spots_per_cell.empty else 0,
                'Median_Spots_Per_Cell': spots_per_cell.median() if not spots_per_cell.empty else 0
            })
            
        except Exception as e:
            print(f"Error processing {base_name}: {str(e)}")
            continue
    
    # Save summary data
    if summary_data:
        pd.DataFrame(summary_data).to_csv(os.path.join(output_dir, 'analysis_summary.csv'), index=False)
        pd.DataFrame(combined_roi_stats).to_csv(os.path.join(output_dir, 'combined_roi_stats.csv'), index=False)

# Set paths for processing
spots_dir = "/Users/ryanpalaganas/Desktop/Areas_of_Responsibility/20-29_Side_Projects/20_AD_TAU_TANGLES/ANALYSIS/max_projections/APLP2_CALM1/split_channels/APLP2/spots_full_field_case_gene/"
masks_dir = "/Users/ryanpalaganas/Desktop/Areas_of_Responsibility/20-29_Side_Projects/20_AD_TAU_TANGLES/ANALYSIS/segmentation_maps/"
output_dir = "/Users/ryanpalaganas/Desktop/Areas_of_Responsibility/20-29_Side_Projects/20_AD_TAU_TANGLES/ANALYSIS/max_projections/APLP2_CALM1/split_channels/APLP2/spots_filtered_case_gene/"

# Run main processing function
process_and_plot_spots(spots_dir, masks_dir, output_dir)
