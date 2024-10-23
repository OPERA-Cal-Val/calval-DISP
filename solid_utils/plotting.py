# Author: Marin Govorcin
# June, 2024
# functions added by Jinwoo, Aug 2024

import pandas as pd
import numpy as np
import matplotlib.patches as patches
from matplotlib import pyplot as plt
from numpy.typing import NDArray
import seaborn as sns

def display_validation(pair_distance: NDArray, pair_difference: NDArray,
                       site_name: str, start_date: str, end_date: str,
                       requirement: float = 2, distance_rqmt: list = [0.1, 50],
                       n_bins: int = 10, threshold: float = 0.683, 
                       sensor:str ='Sentinel-1', validation_type:str='secular',
                       validation_data:str='GNSS'):

   '''
    Parameters:
      pair_distance : array      - 1d array of pair distances used in validation
      pair_difference : array    - 1d array 0f pair double differenced velocity residuals
      site_name : str            - name of the cal/val site
      start_date  : str          - data record start date, eg. 20190101
      end_date : str             - data record end date, eg. 20200101
      requirement : float        - value required for test to pass
                                    e.g, 2 mm/yr for 3 years of data over distance requiremeent
      distance_rqmt : list       - distance over requirement is tested, eg. length scales of 0.1-50 km
      n_bins : int               - number of bins
      threshold : float          - threshold represents percentile of Gaussian normal distribution
                                    within residuals are expected to be to pass the test
                                    e.g. 0.683 for 68.3% or 1-sigma limit 
      sensor : str               - sensor used in validation, e.g Sentinel-1 or NISAR
      validation_type : str      - type of validation: secular, coseismic, transient
      validation_data : str      - data used to validate against; GNSS or INSAR

   Return
      validation_table
      validation_figure
   '''
   # init dataframe
   df = pd.DataFrame(np.vstack([pair_distance,
                                pair_difference]).T,
                                columns=['distance', 'double_diff'])

   # remove nans
   df_nonan = df.dropna(subset=['double_diff'])
   bins = np.linspace(*distance_rqmt, num=n_bins+1)
   bin_centers = (bins[:-1] + bins[1:]) / 2
   binned_df = df_nonan.groupby(pd.cut(df_nonan['distance'], bins),
                                observed=False)[['double_diff']]

   # get binned validation table 
   validation = pd.DataFrame([])
   validation['total_count[#]'] = binned_df.apply(lambda x: np.ma.masked_invalid(x).count())
   validation['passed_req.[#]'] = binned_df.apply(lambda x: np.count_nonzero(x < requirement))
   
   # Add total at the end
   validation = pd.concat([validation, pd.DataFrame(validation.sum(axis=0)).T])
   validation['passed_pc'] = validation['passed_req.[#]'] / validation['total_count[#]']
   validation['success_fail'] = validation['passed_pc'] > threshold
   validation.index.name = 'distance[km]'
   # Rename last row
   validation.rename({validation.iloc[-1].name:'Total'}, inplace=True)

   # Figure
   fig, ax = plt.subplots(1, figsize=(9, 3), layout="none", dpi=200)
   
   # Plot residuals
   ms = 8 if pair_difference.shape[0] < 1e4 else 0.3
   alpha = 0.6 if pair_difference.shape[0] < 1e4 else 0.2
   ax.scatter(df_nonan.distance, df_nonan.double_diff,
              color='black', s=ms, zorder=1, alpha=alpha, edgecolor='None')

   ax.fill_between(distance_rqmt, 0, requirement, color='#e6ffe6', zorder=0, alpha=0.6)
   ax.fill_between(distance_rqmt, requirement, 21, color='#ffe6e6', zorder=0, alpha=0.6)
   ax.vlines(bins, 0, 21, linewidth=0.3, color='gray', zorder=1)
   ax.axhline(requirement, color='k', linestyle='--', zorder=3)

   # Bar plot for each bin
   quantile_th = binned_df.quantile(q=threshold)['double_diff'].values
   for bin_center, quantile, flag in zip(bin_centers,
                                         quantile_th,
                                         validation['success_fail']):
      if flag:
         color = '#227522'
      else:
         color = '#7c1b1b'
      ax.bar(bin_center, quantile, align='center', width=np.diff(bins)[0],
            color='None', edgecolor=color, linewidth=2, zorder=3)
      
   # Add legend with data info
   legend_kwargs = dict(transform=ax.transAxes, verticalalignment='top')
   props = dict(boxstyle='square', facecolor='white', alpha=1, linewidth=0.4)
   textstr = f'Sensor: {sensor} \n{validation_data}-DISP-S1 point pairs\n'
   textstr += f'Record: {start_date}-{end_date}'

   # place a text box in upper left in axes coords
   ax.text(0.02, 0.95, textstr, fontsize=8, bbox=props, **legend_kwargs)
   
   # Add legend with validation info 
   textstr = f'{validation_type.capitalize()} requirement\n'
   textstr += f'Site: {site_name}\n'
   if validation.loc['Total']['success_fail']:
      validation_flag = 'PASSED'
      validation_color = '#239d23'
   else: 
      validation_flag ='FAILED'
      validation_color = '#bc2e2e'

   props = {**props, **{'facecolor':'none', 'edgecolor':'none'}}
   ax.text(0.818, 0.93, textstr, fontsize=8, bbox=props, **legend_kwargs)
   ax.text(0.852, 0.82,  f"{validation_flag}",
           fontsize=10, weight='bold',
           bbox=props, **legend_kwargs)

   rect = patches.Rectangle((0.8, 0.75), 0.19, 0.2,
                           linewidth=1, edgecolor='black',
                           facecolor=validation_color,
                           transform=ax.transAxes)
   ax.add_patch(rect)

   # Title & labels
   fig.suptitle(f"{validation_type.capitalize()} requirement: {site_name}", fontsize=10)
   ax.set_xlabel("Distance (km)", fontsize=8)
   if validation_data == 'GNSS':
       txt = "Double-Differenced \nVelocity Residual (mm/yr)"
   else:
       txt = "Relative Velocity measurement (mm/yr)"    
   ax.set_ylabel(txt, fontsize=8)
   ax.minorticks_on()
   ax.tick_params(axis='x', which='minor', length=4, direction='in', top=False, width=1.5)
   ax.tick_params(axis='both', labelsize=8)
   ax.set_xticks(bin_centers, minor=True)
   ax.set_xticks(np.arange(0,55,5))
   ax.set_ylim(0,20)
   ax.set_xlim(*distance_rqmt)

   validation = validation.rename(columns={'success_fail': f'passed_req [>{threshold*100:.1f}%]'})

   return validation, fig

def display_validation_table(validation_table):
    # Display Statistics
    def bold_last_row(row):
        is_total = row.name == 'Total'
        styles = ['font-weight: bold; font-size: 14px; border-top: 3px solid black' if is_total else '' for _ in row]
        return styles
    
    def style_success_fail(value):
        color = '#e6ffe6' if value else '#ffe6e6'
        return 'background-color: %s' % color

    # Overall pass/fail criterion
    if validation_table.loc['Total'][validation_table.columns[-1]]:
        print("This velocity dataset passes the requirement.")
    else:
        print("This velocity dataset does not pass the requirement.")

    return (validation_table.style
            .bar(subset=['passed_pc'], vmin=0, vmax=1, color='gray')
            .format(lambda x: f'{x*100:.0f}%', na_rep="none", precision=1, subset=['passed_pc'])
            .apply(bold_last_row, axis=1)
            .map(style_success_fail, subset=[validation_table.columns[-1]])
           )

def plot_transient_table(ratio_pd, site, percentage, thresthod, title_text, figname, annot=True):
    _df = ratio_pd.applymap(lambda x: pd.to_numeric(x, errors='coerce'))
    _df.index = _df.index.map(lambda x: f"{x[9:]}")

    # Create a boolean mask for coloring
    # True (1) if value >= 0.5 or not NaN, False (0) otherwise
    _mask = _df.applymap(lambda x: 1 if (x >= thresthod and not pd.isna(x)) else 0)

    # Transpose to swap x and y axes
    _df_transposed = _df.T
    _mask_transposed = _mask.T

    fig, ax = plt.subplots(figsize=(18, 10))

    # Create a custom colormap: red for 0, green for 1
    cmap = plt.cm.colors.ListedColormap(['red', 'green'])

    # Create the heatmap
    if annot:		# heat map has values in cells
        sns.heatmap(_df_transposed, ax=ax, annot=True, cbar=False, cmap=plt.cm.colors.ListedColormap(['white']))
    sns.heatmap(_mask_transposed, cmap=cmap, ax=ax, cbar=False, linewidths=0.5, linecolor='black', alpha=0.2)

    if percentage >= 0.70:
        validation_flag = 'PASSED'
        flag_color = 'green'
    else: 
        validation_flag = 'FAILED'
        flag_color = 'red'
    
    title = ax.set_title(title_text, fontsize=16)

    # Function to create text with mixed font weights
    def mixed_weight_text(text):
        parts = []
        for word in text.split():
            if word == "PASSED" or word == "FAILED":
                parts.append(f"$\\bf{{{word}}}$")  # Bold
            else:
                parts.append(word)
        return " ".join(parts)
    
    # Create the text with mixed weights
    text = mixed_weight_text(f"Transient requirement (Site {site}): {validation_flag}")

    # Get the position of the title
    title_pos = title.get_position()

    # Calculate the position for the text box (just to the right of the title)
    text_x = title_pos[0] + 0.1  # Slightly to the right of the title's center
    text_y = title_pos[1] + 0.03  # Same y-position as the title

    # # Add the text box
    text_box = ax.text(text_x, text_y, text, ha='left', va='center',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor=flag_color, alpha=0.4),
                    transform=ax.transAxes, fontsize=16)

    # Adjust the title's position to make room for the text box
    title.set_position((title_pos[0] - 0.2, title_pos[1]))

    plt.xlabel('Dates', fontsize=16, fontweight ='bold')
    plt.ylabel('Distance Ranges (km)', fontsize=16, fontweight ='bold')

    # Rotate x-axis labels
    plt.xticks(rotation=90)

    ax.axhline(_df_transposed.shape[0]-1, color='black', linewidth=3, linestyle='--')
    if ratio_pd.columns[-1] == 'mean':
        ax.axhline(_df_transposed.shape[0]-2, color='black', linewidth=3, linestyle='--') 

    # Rotate y-axis labels
    plt.yticks(rotation=0)

    # Add a frame around the heatmap
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(2)

    plt.tight_layout()
  
    fig.savefig(figname, bbox_inches='tight', transparent=True, dpi=300)
