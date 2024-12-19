import pandas as pd
import numpy as np
import os
import glob
import tifffile
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries
from scipy.ndimage import center_of_mass

def find_matching_mask(tau_filename, mask_files):
   """
   Find matching segmentation mask for tau image
   """
   parts = tau_filename.split('_')
   ac_index = parts.index('YS') ##change this
   image_number = parts[ac_index - 1]
   orig_part = next(part for part in parts if part.startswith('orig'))
   suffix = orig_part.replace('orig', '')
   
   mask_pattern = f"MAX_{image_number}_YS_NT{suffix}_seg.npy" ##change this
   print(f"Image number: {image_number}")
   print(f"Looking for mask pattern: {mask_pattern}")
   
   matching_masks = [m for m in mask_files if os.path.basename(m) == mask_pattern]
   if matching_masks:
       return matching_masks[0]
   return None

def plot_rois_on_tau(tau_image, mask, output_path, title=""):
    """
    Plot ROIs overlaid on tau image with ROI labels and autoscaled intensity
    """
    plt.figure(figsize=(15, 12))
    
    # Autoscale tau image
    p2, p98 = np.percentile(tau_image, (2, 98))
    tau_scaled = np.clip(tau_image, p2, p98)
    tau_scaled = (tau_scaled - p2) / (p98 - p2)
    
    # Plot scaled tau image in grayscale
    plt.imshow(tau_scaled, cmap='gray')
    
    # Overlay ROI boundaries
    boundaries = find_boundaries(mask)
    plt.imshow(boundaries, cmap='tab20', alpha=0.5)
    
    # Add ROI labels
    unique_rois = np.unique(mask)
    unique_rois = unique_rois[unique_rois > 0]
    
    for roi in unique_rois:
        roi_mask = mask == roi
        center = center_of_mass(roi_mask)
        plt.text(center[1], center[0], str(int(roi)), 
                color='white', fontsize=8, ha='center', va='center',
                bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', pad=1))
    
    plt.title(title)
    plt.axis('off')
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

def analyze_tau_in_rois(tau_image_path, mask_path, output_dir, tau_filename):
   """
   Analyze tau intensity within each ROI and generate overlay plot
   """
   # Read image and mask
   tau_image = tifffile.imread(tau_image_path)
   mask = np.load(mask_path, allow_pickle=True).item()['masks']
   
   # Generate overlay plot
   plot_path = os.path.join(output_dir, tau_filename.replace('.tif', '_overlay.png'))
   plot_rois_on_tau(tau_image, mask, plot_path, title=tau_filename)
   
   # Calculate tau stats for each ROI
   roi_stats = []
   unique_rois = np.unique(mask)
   unique_rois = unique_rois[unique_rois > 0]
   
   for roi in unique_rois:
       roi_mask = mask == roi
       tau_pixels = tau_image[roi_mask]
       
       stats = {
           'ROI_ID': roi,
           'Mean_Tau': np.mean(tau_pixels),
           'Median_Tau': np.median(tau_pixels),
           'Max_Tau': np.max(tau_pixels),
           'Total_Tau': np.sum(tau_pixels),
           'Area_Pixels': np.sum(roi_mask), 
           'Mean_tau_perpixel': (np.mean(tau_pixels) / np.sum(roi_mask)) * 100
       }
       roi_stats.append(stats)
   
   return roi_stats

def process_tau_intensity(tau_dir, masks_dir, output_dir):
   """
   Process all tau images and analyze intensity in ROIs
   """
   os.makedirs(output_dir, exist_ok=True)
   
   tau_files = glob.glob(os.path.join(tau_dir, "*_MAX_C1.tif"))
   mask_files = glob.glob(os.path.join(masks_dir, "*_seg.npy"))
   
   print(f"Found {len(tau_files)} tau images and {len(mask_files)} mask files")
   
   combined_stats = []
   
   for tau_file in tau_files:
       tau_filename = os.path.basename(tau_file)
       print(f"\nProcessing tau image: {tau_filename}")
       
       mask_file = find_matching_mask(tau_filename, mask_files)
       if mask_file is None:
           print(f"No matching mask found for {tau_filename}")
           continue
       
       print(f"Found matching mask: {os.path.basename(mask_file)}")
       
       try:
           parts = tau_filename.split('_')
           ac_index = parts.index('YS') ##change this
           image_number = parts[ac_index - 1]
           orig_part = next(part for part in parts if part.startswith('orig'))
           
           # Updated function call to pass output_dir and filename
           roi_stats = analyze_tau_in_rois(tau_file, mask_file, output_dir, tau_filename)
           
           for stat in roi_stats:
               stat['Image_Number'] = image_number
               stat['Orig_Number'] = f"{image_number}_YS_{orig_part}" ##change this
               stat['ROI_Key'] = f"{image_number}_YS_{orig_part}_ROI{stat['ROI_ID']}" ##change this
           
           combined_stats.extend(roi_stats)
           
       except Exception as e:
           print(f"Error processing {tau_filename}: {str(e)}")
           continue
   
   if combined_stats:
       df = pd.DataFrame(combined_stats)
       output_path = os.path.join(output_dir, 'tau_intensity_analysis.csv')
       df.to_csv(output_path, index=False)
       print(f"\nResults saved to: {output_path}")
       print("\nSummary of tau analysis:")
       print(df.describe())


tau_dir = "/Users/ryanpalaganas/Desktop/Areas_of_Responsibility/20-29_Side_Projects/20_AD_TAU_TANGLES/ANALYSIS/max_projections/YWHAQ_SYT1/split_channels/tau/"
masks_dir = "/Users/ryanpalaganas/Desktop/Areas_of_Responsibility/20-29_Side_Projects/20_AD_TAU_TANGLES/ANALYSIS/segmentation_maps/"
output_dir = "/Users/ryanpalaganas/Desktop/Areas_of_Responsibility/20-29_Side_Projects/20_AD_TAU_TANGLES/ANALYSIS/max_projections/YWHAQ_SYT1/split_channels/tau/intensity/"

process_tau_intensity(tau_dir, masks_dir, output_dir)
