import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import qiime2 as q2
from skbio import OrdinationResults
from joypy import joyplot
from scipy.optimize import curve_fit
import matplotlib.lines as mlines
from sklearn.model_selection import GroupShuffleSplit, GroupKFold, LeaveOneGroupOut
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error, accuracy_score, f1_score


# mapping dictionary for vizualisation labels
name_mapping = {
    "observed_features": "Observed features",
    "shannon_entropy": "Shannon entropy",
    "pielou_evenness": "Pielou evenness",
    "faith_pd": "Faith phylogenetic diversity",
    "age_days": "Age (days)",
    "sex[T.Male]": "Sex (male)",
    "BCQ_Attunement": "Attuned care score (BCQ)",
    "ASQ_Composite": "Infant behavioral \n development (ASQ)",
    "volatility_braycurtis": "Temporal volatility: \n Bray-Curtis dissimilarity",
    "volatility_jaccard": "Temporal volatility: \n Jaccard similarity index",
    "volatility_unweighted_unifrac": "Temporal volatility: \n Unweighted UniFrac distance",
    "volatility_weighted_unifrac": "Temporal volatility: \n Weighted UniFrac distance",
    "time_since_last_bowel_movement_in_h": "Time since last \n bowel movement (h)",
    "prior_time_spent_awake_in_h": "Prior time spent \n awake (h)",
    "time_since_last_feeding_in_h": "Time since last \n feeding (h)",
    "duration_of_last_sleep_in_h": "Duration of last \n sleep (h)",
    "R2_abundance": "Alpha diversity rhythmicity: \n cosine fit score of \n observed features",
    "R2_abundance_evenness": "Alpha diversity rhythmicity: \n cosine fit score of \n Shannon entropy",
    "R2_evenness": "Alpha diversity rhythmicity: \n cosine fit score of \n Pielou evenness",
    "R2_biodiversity": "Alpha diversity rhythmicity: \n cosine fit score of \n Faith's phylogenetic diversity",
    "std_times_between_feedings_in_h": "Std. dev. of time between \n feedings (h)",
    "R2_bifido": "$\it{Bifidobacterium}$  rhythmicity: \n cosine fit score of the $\it{Bifidobacterium}$ genus",
    "R2_veillo": "$\it{Veillonella}$ rhythmicity: \n cosine fit score of the $\it{Veillonella}$ genus",
    "R2_esche": "$\it{Escherichia-Shigella}$ rhythmicity: \n cosine fit score of the $\it{Escherichia-Shigella}$ genus",
    "R2_bacteroides": "$\it{Bacteroides}$ rhythmicity: \n cosine fit score of the $\it{Bacteroides}$ genus",
    "R2_clostridium": "$\it{Clostridium}$ rhythmicity: \n cosine fit score of the $\it{Clostridium}$ genus"}

# labels in  plots with bold formatting
bold_labels = ["Observed features", "Shannon entropy", "Pielou evenness", "Faith phylogenetic diversity",
               "Temporal volatility: \n Bray-Curtis dissimilarity", "Temporal volatility: \n Jaccard similarity index", 
               "Temporal volatility: \n Unweighted UniFrac distance", "Temporal volatility: \n Weighted UniFrac distance",
               "Alpha diversity rhythmicity: \n cosine fit score of \n observed features",
               "Alpha diversity rhythmicity: \n cosine fit score of \n Shannon entropy",
               "Alpha diversity rhythmicity: \n cosine fit score of \n Pielou evenness",
               "Alpha diversity rhythmicity: \n cosine fit score of \n Faith's phylogenetic diversity",
               "$\it{Bifidobacterium}$  rhythmicity: \n cosine fit score of the $\it{Bifidobacterium}$ genus",
               "$\it{Veillonella}$ rhythmicity: \n cosine fit score of the $\it{Veillonella}$ genus",
               "$\it{Escherichia-Shigella}$ rhythmicity: \n cosine fit score of the $\it{Escherichia-Shigella}$ genus",
               "$\it{Bacteroides}$ rhythmicity: \n cosine fit score of the $\it{Bacteroides}$ genus",
               "$\it{Clostridium}$ rhythmicity: \n cosine fit score of the $\it{Clostridium}$ genus"]


def process_qza_file(qza_file_path, stat_function):
    """
    Processes a .qza file to load and transform data into a processed dataframe.

    Parameters:
    qza_file_path (str): Path to the input .qza file.
    stat_function (function): A function that takes a dataframe as input
                                     and returns a pocessed dataframe.

    Returns:
    pd.DataFrame: The processed dataframe.
    """
    # Load the Qiime 2 Artifact as a dataframe
    final_feature_table_df = load_qza_artifact(qza_file_path)

    # Perform function
    final_feature_table_percent = stat_function(final_feature_table_df)

    return final_feature_table_percent


def compute_stats_for_percentage_columns(df):
    """
    Compute the mean, median, and standard deviation of all columns 
    with names including the string "percentage".
    
    :param df: Input pandas DataFrame
    :return: DataFrame with computed statistics (mean, median, std)
    """
    # Filter columns that include the string "percentage" in their names
    percentage_columns = [col for col in df.columns if "percentage" in col]
    
    # Select only the percentage columns
    df_percentage = df[percentage_columns]
    
    # Compute statistics
    stats = {
        "mean": df_percentage.mean(),
        "median": df_percentage.median(),
        "std": df_percentage.std()
    }
    
    # Convert the dictionary of statistics into a DataFrame
    stats_df = pd.DataFrame(stats)
    
    return stats_df


def load_qza_artifact(path, case='metadata'):
    """
    Loads a .qza artifact and converts it into a pandas DataFrame.

    Parameters:
    path (str): Path to the .qza file.
    case (str): Type of artifact to load ('metadata' or 'ordination').

    Returns:
    pd.DataFrame: Extracted data as a DataFrame.
    """
    artifact = q2.Artifact.load(path)

    if case == 'metadata':
        df = artifact.view(q2.Metadata).to_dataframe()
    elif case == 'ordination':
        df = artifact.view(OrdinationResults)
    else:
        raise ValueError(f"Unsupported case: {case}. Choose 'metadata' or 'ordination'.")

    return df


def process_pcoa(pcs, metadata):
    """
    Extracts the first three PCoA axes and merges with metadata.

    Parameters:
    pcs (pd.DataFrame): The PCoA results containing samples.
    metadata (pd.DataFrame): The metadata to merge.

    Returns:
    pd.DataFrame: Processed DataFrame with the first 3 PCoA axes and metadata.
    """
    pcs_processed = pcs.samples.iloc[:, :3]  # Take first 3 columns
    pcs_processed.columns = ['Axis 1', 'Axis 2', 'Axis 3']  # Rename columns
    return pd.concat([pcs_processed, metadata], join="inner", axis=1)  # Merge with metadata


def load_tsv(file_path, index_col=0):
    """
    Loads a .tsv file into a Pandas DataFrame.

    Parameters:
    file_path (str): Path to the .tsv file.
    index_col (int or None): Column to use as the row labels. Default is 0 (first column).

    Returns:
    pd.DataFrame: The loaded DataFrame.
    """
    return pd.read_csv(file_path, sep='\t', index_col=index_col)


def save_df_as_tsv(df, file_path, index_name="id"):
    """
    Renames the index, saves a DataFrame as a .tsv file, and ensures proper formatting.

    Parameters:
    df (pd.DataFrame): The DataFrame to be saved.
    file_path (str): The output file path, including the filename (e.g., 'output/data.tsv').
    index_name (str): The name to assign to the index column (default: 'id').

    Returns:
    None
    """
    df = df.rename_axis(index_name)  # Set index name
    df.to_csv(file_path, sep='\t', index=True)  # Save as .tsv


def get_taxonomy(col_name, level='p'):
    """
    Extracts a specific taxonomic level from a taxonomy string.

    Parameters:
    col_name (str): The taxonomy string (e.g., 'k__Bacteria; p__Firmicutes').
    level (str): The prefix to search for (e.g., 'p' for phylum, 'g' for genus).

    Returns:
    str: The identified taxonomic rank or a fallback 'Unassigned' string.
    """
    # Map prefixes to names for the return message
    ranks = {'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus'}
    target = f"{level}__"
    
    # Split the taxonomy string by semicolon
    parts = col_name.split(';')

    # Find the part that starts with the target
    for part in parts:
        clean_part = part.strip()
        if clean_part.startswith(target):
            return clean_part
            
    return f"Unassigned {ranks.get(level, 'rank')}"


def get_tax_per_timepoint(tax_ft, metdata_ft, timepoint_col):
    """
    Helper function to get the relative abundance of a taxonomic level per timepoint.
    """

    # Merge the taxonomic feature table with metadata
    df_with_metadata = tax_ft.merge(metdata_ft[timepoint_col], left_index=True, right_index=True)
    # Group by timepoint and sum
    df_with_metadata_timepoint = df_with_metadata.groupby(timepoint_col).sum()
    # Convert to relative abundance
    df_with_metadata_timepoint = df_with_metadata_timepoint.div(df_with_metadata_timepoint.sum(axis=1), axis=0)
    # Check that rows sum to 1
    row_sums = df_with_metadata_timepoint.sum(axis=1)
    if not np.all(np.isclose(row_sums, 1.0)):
        print("Warning: Some total rel. abund. do not sum to 1")

    return df_with_metadata_timepoint


def collapse_taxa(df, min_rel_abundance=0, top_n=10, keep_col='Unassigned'):
    """
    Collapses a taxonomy DataFrame by retaining top taxa and grouping the rest into 'Others'.

    Parameters:
    df (pd.DataFrame): The DataFrame to be collapsed (samples as rows, taxa as columns).
    min_rel_abundance (float): The minimum abundance threshold for a taxon to be kept.
    top_n (int): The number of most abundant taxa to keep if min_rel_abundance is 0.
    keep_col (str): The name of a specific column to preserve from the ranking logic (default: 'Unassigned').

    Returns:
    pd.DataFrame: A DataFrame containing the top taxa, the preserved column, and a new 'Others' column.
    """
    # 1. Identify columns to rank (exclude 'Unassigned' from the ranking logic)
    cols_to_rank = [c for c in df.columns if c != keep_col]
    
    if min_rel_abundance > 0:
        # Filter columns based on the minimum relative abundance threshold
        print("Filtering taxa based on the minimum relative abundance threshold: ", min_rel_abundance)
        top_taxa = df.columns[df.max() > min_rel_abundance]
    else:    
        # 2. Get the top N taxa based on the sum across all rows
        print("Selecting top N taxa based on total abundance: ", top_n)
        top_taxa = df[cols_to_rank].sum().sort_values(ascending=False).head(top_n).index.tolist()
    
    # 3. Identify columns that are not in the top N and not the protected 'Unassigned' column
    others_cols = [c for c in df.columns if c not in top_taxa and c != keep_col]
    
    # 4. Build the new DataFrame
    # Start with the top taxa
    top_taxa_sorted = df[top_taxa].sum().sort_values(ascending=False).index.tolist()
    df_collapsed = df[top_taxa_sorted].copy()
    
    # Add 'Unassigned' if it exists in the original data
    if keep_col in df.columns:
        df_collapsed[keep_col] = df[keep_col]
    
    # Sum the remaining columns into 'Others'
    df_collapsed['Others'] = df[others_cols].sum(axis=1)

    df_final = df_collapsed.copy()
    
    return df_final


def plot_taxa_stacked_bar(df, level_name, out_path, threshold=0, top_n=10, figsize=(5.5, 4)):
    """
    Plots a stacked bar chart of relative abundance.
    
    Parameters:
    - df: The prepared DataFrame (rows=timepoints, cols=taxa)
    - level_name: String (e.g., 'Phylum' or 'Genus') for labels/titles
    - out_path: Directory path to save the PDF
    - threshold: Minimum max-abundance to include a taxon (default 1%)
    - top_n: Number of top taxa to include
    - figsize: Figure size tuple
    """
    
    if level_name == "phylum":
        keep_col = 'Unassigned phylum'
    elif level_name == "genus":
        keep_col = 'Unassigned genus'
    else:
        keep_col = 'Unassigned'

    # Filter for abundance threshold
    plot_df = collapse_taxa(df, min_rel_abundance=threshold, top_n=top_n, keep_col=keep_col)
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    # plotting
    ax.stackplot(
                plot_df.index,
                plot_df.T,
                labels=plot_df.columns,
                colors=plt.cm.tab20.colors)

    # Legend handling: reverse order to match the visual stack (top-to-bottom)
    handles, labels = ax.get_legend_handles_labels()
    labels = [label.split('__')[-1] if '__' in label else label for label in labels]

    ax.legend(handles[::-1], labels[::-1], title=level_name, 
              loc='upper left', bbox_to_anchor=(1, 1))

    # Aesthetics
    ax.tick_params(which="both", direction="out")
    ax.set_ylabel("Relative abundance")
    ax.set_title(f"Relative abundance of gut bacterial {level_name} over time")

    fig.tight_layout()
    
    # Save and show
    save_filename = f"{out_path}/{level_name.lower()}_relative_abundance_over_time.pdf"
    fig.savefig(save_filename, dpi=300, bbox_inches='tight')
    plt.show()


def compute_distance_to_centroid(df, new_col_name, x_col='Axis 1', y_col='Axis 2', z_col='Axis 3', 
                                 centroid_x_col='Centroid x', centroid_y_col='Centroid y', centroid_z_col='Centroid z'):
    """
    Computes the Euclidean distance from each point (x, y, z) to its centroid (centroid_x, centroid_y, centroid_z).
    
    Parameters:
    df (pd.DataFrame): DataFrame containing point and centroid coordinates.
    x_col, y_col, z_col (str): Column names for the point coordinates.
    centroid_x_col, centroid_y_col, centroid_z_col (str): Column names for the centroid coordinates.
    
    Returns:
    pd.DataFrame: Original DataFrame with an additional column 'distance_to_centroid'.
    """
    df[new_col_name] = np.sqrt(
        (df[x_col] - df[centroid_x_col])**2 +
        (df[y_col] - df[centroid_y_col])**2 +
        (df[z_col] - df[centroid_z_col])**2
    )
    return df[new_col_name]

    
def min_max_normalize(series, min_val, max_val, reverse=False):
    """Normalize a value using min-max scaling to range [0,1]. Optionally reverse the scale."""
    normalized = (series - min_val) / (max_val - min_val)
    return 1 - normalized if reverse else normalized

    
def compute_sleep_score(dataframe, weights=None):
    """
    Compute a composite sleep score for each row in a dataframe.
    
    Parameters:
    - df: DataFrame with columns ['nighttime_sleep_duration_in_h', 'sleep_latency_in_h', 'bedtime_in_h', 'number_of_awakenings']
    - weights: Dictionary with weights for each variable (default: equal weights)
    
    Returns:
    - Additional column 'sleep_score'
    """
    df = dataframe.copy()

    sleep_dur = min_max_normalize(df['nighttime_sleep_duration_in_h'], 
                                  df['nighttime_sleep_duration_in_h'].min(), 
                                  df['nighttime_sleep_duration_in_h'].max())
    sleep_lat = min_max_normalize(df['sleep_latency_in_h'], 
                                  df['sleep_latency_in_h'].min(), 
                                  df['sleep_latency_in_h'].max(), reverse=True)
    bedtime = min_max_normalize(df['bedtime_in_h'], 
                                df['bedtime_in_h'].min(), 
                                df['bedtime_in_h'].max(), reverse=True)
    awakenings = min_max_normalize(df['number_of_awakenings'], 
                                   df['number_of_awakenings'].min(), 
                                   df['number_of_awakenings'].max(), reverse=True)
    
    if weights is None:
        weights = { 'nighttime_sleep_duration_in_h': 0.25, 
                   'sleep_latency_in_h': 0.25, 
                   'bedtime_in_h': 0.25, 
                   'number_of_awakenings': 0.25 }
    
    # Compute weighted sum
    df['composite'] = (weights['nighttime_sleep_duration_in_h'] * sleep_dur +
                            weights['sleep_latency_in_h'] * sleep_lat +
                            weights['bedtime_in_h'] * bedtime +
                            weights['number_of_awakenings'] * awakenings)
    
    return df['composite']


def cosine_function(x, A, phi, C):
    """
    Cosine function to model rhythmicity. Here lambda/period = 24h.
    """
    return A * np.cos((2 * np.pi / 24) * x + phi) + C


def fit_cosine_and_plot(data, timepoint_values, timepoints_colors, feature_column, ax, show_legend=True):
    """
    Fits a cosine function to the given feature column for each timepoint and subject.
    
    Parameters:
    - data: DataFrame containing the dataset.
    - timepoint_values: List of timepoints to analyze.
    - timepoints_colors: List of colors corresponding to timepoints.
    - feature_column: Column name to analyze (e.g., "observed_features", "shannon_entropy").
    - ax: Matplotlib Axes object to plot on.
    - show_legend: Boolean, whether to display the legend.

    Returns:
    - fit_scores_df: DataFrame containing R², RMSE, and MAE for each subject.
    - summary_stats: Summary statistics (mean, median) per timepoint.
    """

    # Store results
    subject_fit_scores = []
    timepoint_counts = {}  # Store (infants, samples) per timepoint
    legend_handles = {}

    for i, timepoint in enumerate(timepoint_values):
        subset_data_time = data[data['timepoint'] == timepoint]
        color = timepoints_colors[i]

        n_infants = 0
        n_samples = 0

        for subject in subset_data_time['infant_id'].unique():
            subset_data = subset_data_time[subset_data_time['infant_id'] == subject]

            if subset_data.empty or len(subset_data) < 4:
                continue  # Skip subjects with too few points

            n_infants += 1
            n_samples += len(subset_data)

            X = np.array(subset_data['hour_stool_sample'])
            y = np.array(subset_data[feature_column])

            # Initial parameter guesses
            # Amplitude (A): Half the range of y_data.
            # Phase (φ): Start at 0.
            # Offset (C): Mean of y_data.
            A_init = (np.max(y) - np.min(y)) / 2
            phi_init = 0
            C_init = np.mean(y)
            initial_guess = [A_init, phi_init, C_init]

            params, _ = curve_fit(cosine_function, X, y, p0=initial_guess, maxfev=5000,
                                  bounds=([0, -np.pi, np.min(y)], [np.inf, np.pi, np.max(y)]))
            A_fit, phi_fit, C_fit = params

            # Generate fit curve
            X_fit = np.linspace(0, 24, 100)
            y_fit = cosine_function(X_fit, A_fit, phi_fit, C_fit)

            scatter_handle = ax.scatter(X, y, color=color, alpha=0.6)
            ax.plot(X_fit, y_fit, color=color, linewidth=2)

            if timepoint not in legend_handles:
                legend_handles[timepoint] = scatter_handle

            # Compute fit quality metrics
            r2 = r2_score(y, cosine_function(X, A_fit, phi_fit, C_fit))
            rmse = np.sqrt(mean_squared_error(y, cosine_function(X, A_fit, phi_fit, C_fit)))
            mae = mean_absolute_error(y, cosine_function(X, A_fit, phi_fit, C_fit))

            # Store results
            subject_fit_scores.append({"infant_id": subject,
                                        "timepoint": timepoint,
                                        "R²": r2,
                                        "RMSE": rmse,
                                        "MAE": mae,
                                        "amplitude": A_fit,
                                        "constant": C_fit})

        timepoint_counts[timepoint] = (n_infants, n_samples)

    # Convert results to DataFrame
    fit_scores_df = pd.DataFrame(subject_fit_scores)

    # Compute summary stats per timepoint
    summary_stats = fit_scores_df.groupby("timepoint").agg({
        "R²": ["mean", "median"],
        "RMSE": ["mean", "median"],
        "MAE": ["mean", "median"],
        "amplitude": "median",
        "constant": "median"
    }).reset_index()

    legend_labels = [f"{tp} ({timepoint_counts[tp][0]} infants, {timepoint_counts[tp][1]} samples)" 
                     for tp in legend_handles.keys()]

    ax.set_title(f"Cosine fit of {feature_column.replace('_', ' ')}", fontsize=16)
    ax.set_xlabel("hour_stool_sample", fontsize=14)
    ax.set_ylabel(feature_column.replace("_", " ").title(), fontsize=14)
    ax.grid(True)
    ax.set_xticks(range(0, 25, 4))

    if show_legend:
        ax.legend(legend_handles.values(), legend_labels, title="", bbox_to_anchor=(1, 1), loc='upper left', fontsize=14)

    return fit_scores_df, summary_stats


def fit_cosine_and_plot_scaled(data, timepoint_values, timepoints_colors, feature_column, 
                               ax, fit_scores, show_legend=True, median_amplitudes=None, median_constants=None):
    """
    Plots a scaled cosine fit using precomputed median amplitude and constant per timepoint.

    Parameters:
    - data: DataFrame containing the dataset.
    - timepoint_values: List of unique timepoints.
    - timepoints_colors: List of colors corresponding to timepoints.
    - feature_column: Column to analyze (e.g., "observed_features", "shannon_entropy").
    - ax: Matplotlib Axes object to plot on.
    - show_legend: Boolean, whether to display the legend.
    - median_amplitudes: Precomputed dictionary of median amplitude per timepoint.
    - median_constants: Precomputed dictionary of median constant per timepoint.

    Returns:
    - None
    """

    if median_amplitudes is None or median_constants is None:
        raise ValueError("median_amplitudes and median_constants must be provided.")
    
    timepoint_counts = {}  # Store (infants, samples) per timepoint
    legend_handles = {}

    # Plot scaled version using precomputed median values
    for i, timepoint in enumerate(timepoint_values):
        subset_data_time = data[data['timepoint'] == timepoint]
        subset_R2_time = fit_scores[fit_scores['timepoint'] == timepoint]
        color = timepoints_colors[i]
        n_infants = 0
        n_samples = 0

        for subject in subset_data_time['infant_id'].unique():
            subset_data = subset_data_time[subset_data_time['infant_id'] == subject]
            subset_R2 = subset_R2_time[subset_R2_time['infant_id'] == subject]

            if subset_data.empty or len(subset_data) < 4:
                continue
            
            if subset_R2["R²"].iloc[0] < 0.5:
                continue

            n_infants += 1
            n_samples += len(subset_data)

            X = np.array(subset_data['hour_stool_sample'])
            y = np.array(subset_data[feature_column])

            # Initial parameter guesses:
            # Amplitude (A): Half the range of y_data.
            # Phase (φ): Start at 0.
            # Offset (C): Mean of y_data.
            A_init = (np.max(y) - np.min(y)) / 2
            phi_init = 0
            C_init = np.mean(y)
            initial_guess = [A_init, phi_init, C_init]

            params, _ = curve_fit(cosine_function, X, y, p0=initial_guess, maxfev=5000,
                              bounds=([0, -np.pi, np.min(y)], [np.inf, np.pi, np.max(y)]))
            A_fit, phi_fit, C_fit = params

            X_fit = np.linspace(0, 24, 100)
            A_fit = list(median_amplitudes)[i][timepoint]
            C_fit = list(median_constants)[i][timepoint]
            y_fit = cosine_function(X_fit, A_fit, phi_fit, C_fit)

            scatter_handle = ax.scatter(X, y, color=color, alpha=0.6)
            ax.plot(X_fit, y_fit, color=color, linewidth=2)

            if timepoint not in legend_handles:
                legend_handles[timepoint] = scatter_handle
    
        timepoint_counts[timepoint] = (n_infants, n_samples)
    
    legend_labels = [f"{tp} ({timepoint_counts[tp][0]} infants, {timepoint_counts[tp][1]} samples)" for tp in legend_handles.keys()]

    ax.set_title(f"Cosine fit (R² > 0.5) of {feature_column.replace('_', ' ')} \n - Adjusted by Median Amplitude & Constant", fontsize=16)
    ax.set_xlabel("hour_stool_sample", fontsize=14)
    ax.set_ylabel(feature_column.replace("_", " ").title(), fontsize=14)
    ax.grid(True)
    ax.set_xticks(range(0, 25, 4))

    if show_legend:
        ax.legend(legend_handles.values(), legend_labels, title="", loc='lower center', fontsize=14)

        
def blend_with_white(color, blend_factor=0.3):
    """
    Blends a given RGB color with white.
    `blend_factor=0.3` means 30% white, 70% original color.
    """
    white = np.array([1, 1, 1, 1])  # RGBA white
    return tuple((1 - blend_factor) * np.array(color) + blend_factor * white)


def samples_scatter_plot(data_analysis, data_all, timepoints, output_file):
    """
    Creates a scatter plot of sample collection times with vertical bars for additional samples.

    Parameters:
        data_analysis (DataFrame): Data containing 'hour_stool_sample', 'timepoint', and 'infant_id'.
        data_all (DataFrame): Data containing additional samples with 'hour_stool_sample', 'timepoint', and 'infant_id'.
        timepoints (list): List of timepoints to include in the plot (e.g., ['2 months', '4 months', '6 months']).
        output_file (str): Name of the file to save the plot.
    """
    num_plots = len(timepoints)

    fig, axes = plt.subplots(1, num_plots, figsize=(6 * num_plots, 6), sharey=True)

    # Sort infants by infant_id
    infants = sorted(data_analysis['infant_id'].unique())  # Sorting the infant IDs
    infants_num = list(range(1, len(data_analysis['infant_id'].unique()) + 1))
    y_positions = range(len(infants))

    for i, timepoint in enumerate(timepoints):
        ax = axes[i]
        timepoint_data_analysis = data_analysis[data_analysis['timepoint'] == timepoint]
        
        # Count valid samples for the current timepoint
        num_samples = timepoint_data_analysis.dropna(subset=['hour_stool_sample']).shape[0]
        
        # Add vertical bars for the samples from `data_all` first, to be behind the scatter points
        timepoint_data_all = data_all[data_all['timepoint'] == timepoint]
        for j, infant in enumerate(infants):
            infant_data_all = timepoint_data_all[timepoint_data_all['infant_id'] == infant]
            valid_data_all = infant_data_all.dropna(subset=['hour_stool_sample'])
            
            # Increase the length of the vertical bars (set y-axis range from j-0.5 to j+0.5)
            ax.vlines(valid_data_all['hour_stool_sample'], j - 0.5, j + 0.5, color='black', alpha=0.4)
        
        # Now plot the scatter plot for the data_analysis (dots will be in front)
        for j, infant in enumerate(infants):
            infant_data_analysis = timepoint_data_analysis[timepoint_data_analysis['infant_id'] == infant]
            valid_data_analysis = infant_data_analysis.dropna(subset=['hour_stool_sample'])
            ax.scatter(valid_data_analysis['hour_stool_sample'], [j] * len(valid_data_analysis), alpha=1)

        ax.set_xticks(range(0, 24, 4))  # Add x-axis ticks for time of day
        ax.set_xlim(0, 24)  # Set x-axis limit from 0 to 24
        ax.set_xlabel('Time of day (hours)', fontsize=14)
        ax.set_title(f'Sample collection times at {timepoint} ({num_samples} samples)', fontsize=14)
        ax.grid(axis='y')

    # Set common y-axis label
    axes[0].set_yticks(y_positions)
    axes[0].set_yticklabels(infants_num)
    axes[0].set_ylabel('Infants', fontsize=14)

    plt.tight_layout()

    # Save and show the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()


def density_plot(data_analysis, timepoints, colors, output_file, alpha=1, ticklabel=14, labelsize=16, titlesize=16):
    """
    Creates a density plot (joyplot) for sample collection times with sample counts included in the labels.

    Parameters:
        data_analysis (DataFrame): Data containing 'hour_stool_sample' and 'timepoint' columns.
        timepoints (list): List of timepoints to include in the plot (e.g., ['2 months', '4 months', '6 months']).
        colors (list): List of colors for the density plots, one per timepoint.
        alpha (float): Transparency level for the density plots.
        ticklabel (int): Font size for tick labels.
        labelsize (int): Font size for axis labels.
        titlesize (int): Font size for the plot title.
        output_file (str): Name of the file to save the plot.
    """
    # Creating a dictionary of data for each timepoint and calculating sample sizes
    timepoint_data = {}
    updated_timepoints = []
    for timepoint in timepoints:
        subset = data_analysis[data_analysis['timepoint'] == timepoint]['hour_stool_sample'].dropna()
        timepoint_data[timepoint] = subset
        n_samples = len(subset)
        updated_timepoints.append(f"(n={n_samples} samples)")

    ls_colors = colors  # List of colors for each timepoint

    # Generate joyplots
    fig, ax = joyplot(
        data=list(timepoint_data.values()),  # The list of sample hour data for each timepoint
        labels=updated_timepoints,  # Updated labels with sample sizes
        overlap=0,  # Control overlap between plots
        alpha=alpha * 0.9,
        ylabelsize=ticklabel,
        xlabelsize=ticklabel,
        color=ls_colors,  # Color for each plot
        linewidth=0.8,  # Line width for the outlines
        figsize=(6, 7),  # Size of the plot
    )

    # Set x-axis limits and ticks
    for axis in ax:
        axis.set_xlim(0, 24)
        axis.set_xticks([0, 6, 12, 18, 24])

    # Set y-axis label
    # Customize plot appearance
    plt.rc("font", size=labelsize)
    plt.title("", fontsize=titlesize, y=0.9)
    plt.xlabel("Time of the day (hours)", fontsize=labelsize)
    plt.grid(True, axis='x')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()


def generate_color_palette(list):
    """
    Generates a custom color palette composed of lightened shades derived from
    selected colors within the Spectral colormap.

    Parameters:
    list (list): A list of three integers specifying the number of shades to generate
                 from each of the three selected base colors.
    """

    # Get colors from the Spectral colormap
    colors = [plt.cm.Spectral(i / float(6)) for i in range(7)]
    timepoints_colors = [colors[0], colors[5], colors[6]]

    # Function to generate shades of a color
    def generate_shades(base_color, num_shades, lighten=True):
        base_rgb = np.array(base_color[:3])  # Ignore alpha channel
        white_or_black = np.array([1, 1, 1]) if lighten else np.array([0, 0, 0])
        
        shades = [tuple((1 - t) * base_rgb + t * white_or_black) for t in np.linspace(0, 1, num_shades)]
        return shades

    # Generate shades
    shades_0 = generate_shades(timepoints_colors[0], list[0], lighten=True)  # Lightening colors[0]
    shades_5 = generate_shades(timepoints_colors[1], list[1], lighten=True)  # Lightening colors[5]
    shades_6 = generate_shades(timepoints_colors[2], list[2], lighten=True)  # Lightening colors[6]

    # Combine all shades
    final_palette = shades_0 + shades_5 + shades_6
    
    return final_palette


def plot_pcoa_infants(data, column, title_suffix, ax, var_explained, final_palette="husl"):
    """
    Function to plot a PCoA scatterplot grouped by infant_id with no specific hue order
    and without displaying the legend.
    
    Parameters:
    - data: DataFrame containing the PCoA data with 'Axis 1', 'Axis 2', and the grouping column.
    - column: Column to use for hue (e.g., 'infant_id').
    - title_suffix: String to be added to the plot title.
    - ax: Axis object on which to plot the graph.
    - var_explained: List containing the percentage of variance explained by Axis 1 and Axis 2.
    - final_palette: Color palette to use for the plot (default is "husl").
    """

    # Get unique colors from the scatterplot

    data = data.sort_values(column)

    unique_infants = data[column].unique()

    if final_palette=="husl":
        palette = sns.color_palette("husl", len(unique_infants))
    else:
        palette = final_palette

    # Create the scatterplot on the given axis
    scatter = sns.scatterplot(data=data,
                               x='Axis 1',
                               y='Axis 2',
                               hue=column,
                               palette=palette,  # Use a diverse color palette
                               ax=ax,
                               legend=False)  # Disable the legend

    # Calculate centroids for each group (infant_id in this case)
    centroids = data.groupby(column)[['Axis 1', 'Axis 2']].mean().reset_index()

    color_map = dict(zip(unique_infants, palette))

    # Overlay filled, round centroids for each infant_id
    for _, row in centroids.iterrows():
        ax.scatter(row['Axis 1'], 
                   row['Axis 2'], 
                   color=color_map[row[column]], 
                   s=300, 
                   alpha=0.6,  # Adjust alpha for shading
                   edgecolor='black',  # Add edge for better visibility
                   linewidth=1.2)

    ax.set_xlabel(f"Axis 1 ({var_explained[0]})")
    ax.set_ylabel(f"Axis 2 ({var_explained[1]})")

    # Add title and labels
    ax.set_title(f"PCoA with distance matrix as {title_suffix}", fontsize=15)
    ax.grid(True)


def violinplot_with_lines(data, x_column, y_column, subject_column, title, y_column_name, palette, 
                          ax=None, fontsize_title=14, save_path=None, text=None, loc=0.05):
    """
    Creates a violin plot with overlaid boxplots and subject-level connecting lines
    across timepoints (or any ordered x-axis variable).

    Parameters
    ----------
    data : pandas.DataFrame
        The dataset containing all required columns.
    x_column : str
        Column name for the categorical x-axis variable (e.g., timepoint).
    y_column : str
        Column name for the numerical y-axis variable.
    subject_column : str
        Column name identifying subjects for drawing repeated-measure lines.
    title : str
        Title of the plot.
    y_column_name : str
        Y-axis label.
    palette : list or dict
        Color palette to assign to each x-axis category for violin and boxplots.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on. If None, a new figure and axis are created.
    fontsize_title : int, optional (default = 14)
        Font size for the plot title. Reduced automatically when using subplots.
    save_path : str, optional
        If provided, saves the figure to this path (only if a new figure is created).
    text : str, optional
        Additional annotation text to display in the top-left corner of the plot.
    loc : float, optional (default = 0.05)
        Horizontal position (axis fraction) of the annotation text.
    """

    # Sort data by x-axis values
    data = data.sort_values([x_column])
    timepoint_order = data[x_column].unique()

    # Get count of each timepoint
    count_infants = data[x_column].dropna().value_counts().sort_index()
    x_labels = [f"{tp} (n={count_infants[tp]})" for tp in count_infants.index]

    # If ax is not provided, create a new figure and axis
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
        created_fig = True  # Track if we created a new figure
    else:
        fontsize_title = 16  # Reduce title size for subplots

    # Create violin plot
    sns.violinplot(x=data[x_column], y=data[y_column], palette=palette, ax=ax, showfliers=False, cut=0)

    # Overlay boxplot
    sns.boxplot(x=data[x_column], y=data[y_column], saturation=0.3, width=0.4, palette=palette, 
                boxprops={'zorder': 2}, ax=ax)

    # Overlay individual subject lines
    x_positions = {tp: i for i, tp in enumerate(timepoint_order)}
    for subject_id, sub_df in data.groupby(subject_column):
        sub_df = sub_df.sort_values(x_column)
        ax.plot([x_positions[tp] for tp in sub_df[x_column]], sub_df[y_column], linestyle='-', 
                alpha=0.2, color='black')

    # Set titles and labels
    ax.set_title(title, fontsize=fontsize_title)
    ax.set_xlabel(None)
    ax.set_ylabel(y_column_name, fontsize=14)
    ax.grid(True)

    # Rotate x-axis labels if using subplots
    rotation_angle = 45 if not created_fig else 0
    ax.set_xticklabels(x_labels, fontsize=12, rotation=rotation_angle, ha="right" if rotation_angle else "center")

    # Add optional text annotation
    if text is not None:
        ax.text(loc, 0.95, text, transform=ax.transAxes, fontsize=9, verticalalignment='top', 
                horizontalalignment='left', bbox=dict(facecolor='white', alpha=0.5))

    # Save the figure if path is provided
    if save_path and created_fig:  
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    # If a new figure was created, show it
    if created_fig:
        plt.show()

    return ax  # Return axis for further modifications if needed


def violinplot_no_lines(data, x_column, y_column, title, y_column_name, palette, 
                          ax=None, fontsize_title=16, fontsize_axes=14, save_path=None, text=None, loc=0.05):
    """
    Creates a violin plot with overlaid boxplots.

    Parameters
    ----------
    data : pandas.DataFrame
        The dataset containing all required columns.
    x_column : str
        Column name for the categorical x-axis variable (e.g., timepoint).
    y_column : str
        Column name for the numerical y-axis variable.
    title : str
        Title of the plot.
    y_column_name : str
        Y-axis label.
    palette : list or dict
        Color palette to assign to each x-axis category for violin and boxplots.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on. If None, a new figure and axis are created.
    fontsize_title : int, optional (default = 16)
        Font size for the plot title. Reduced automatically when using subplots.
    fontsize_axes : int, optional (default = 14)
        Font size for axis labels and tick labels.
    save_path : str, optional
        If provided, saves the figure to this path (only if a new figure is created).
    text : str, optional
        Additional annotation text to display in the top-left corner of the plot.
    loc : float, optional (default = 0.05)
        Horizontal position (axis fraction) of the annotation text.
    """
    
    # Sort data by x-axis values
    data = data.sort_values([x_column])

    # Get count of each timepoint
    count_samples = data[x_column].dropna().value_counts().sort_index()
    x_labels = [f"{tp} (n={count_samples[tp]})" for tp in count_samples.index]

    # If ax is not provided, create a new figure and axis
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
        created_fig = True  # Track if we created a new figure
    else:
        fontsize_title = 16  # Reduce title size for subplots

    # Create violin plot
    sns.violinplot(x=data[x_column], y=data[y_column], palette=palette, ax=ax, cut=0)

    # Overlay boxplot
    sns.boxplot(x=data[x_column], y=data[y_column], saturation=0.3, width=0.4, palette=palette, 
                boxprops={'zorder': 2}, ax=ax)

    # Set titles and labels
    ax.set_title(title, fontsize=fontsize_title)
    ax.set_xlabel(None)
    ax.set_ylabel(y_column_name, fontsize=fontsize_axes)
    ax.grid(True)

    # Rotate x-axis labels if using subplots
    rotation_angle = 45 if not created_fig else 0
    ax.set_xticklabels(x_labels, fontsize=fontsize_axes, rotation=rotation_angle, ha="right" if rotation_angle else "center")

    # Add optional text annotation
    if text is not None:
        ax.text(loc, 0.95, text, transform=ax.transAxes, fontsize=10, verticalalignment='top', 
                horizontalalignment='left', bbox=dict(facecolor='white', alpha=0.5))

    # Save the figure if path is provided
    if save_path and created_fig:  
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    # If a new figure was created, show it
    if created_fig:
        plt.show()

    return ax  # Return axis for further modifications if needed



def scatterplot(data, x_column, y_column, timepoint_column, title, y_column_name, x_column_name, palette, 
                            ax=None, fontsize_title=16, fontsize_axes=14, save_path=None, 
                            loc_position = 'lower right', log=False):
    """
    Creates a scatter plot with regression lines and confidence intervals for different subjects or timepoints. Also annotates the legend with sample sizes.
    
    Parameters
    ----------
    data : pandas.DataFrame
        The dataset containing values for the x-axis, y-axis, and groups defined by `timepoint_column`.
    x_column : str
        Name of the column in `data` to be plotted on the x-axis.
    y_column : str
        Name of the column in `data` to be plotted on the y-axis.
    timepoint_column : str
        Name of the column that defines groups (e.g., timepoints) used for coloring and legend entries.
    title : str
        Title of the plot.
    y_column_name : str
        Display name (label) for the y-axis.
    x_column_name : str
        Display name (label) for the x-axis.
    palette : list or dict
        Color palette used for different timepoints. Should have at least as many entries as unique groups in `timepoint_column`.
    ax : matplotlib.axes.Axes, optional
        A Matplotlib Axes object to draw the plot on. If None, a new figure and axes are created.
    fontsize_title : int, optional (default=16)
        Font size for the plot title.
    fontsize_axes : int, optional (default=14)
        Font size for the axis labels.
    save_path : str or None, optional
        If provided, the figure will be saved to this file path.
    loc_position : str or None, optional (default='lower right')
        Position of the legend. If None, the legend is omitted.
    """

    data = data.sort_values([timepoint_column])

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 6))
    
    # Scatter plot
    sns.scatterplot(data=data, x=x_column, y=y_column, hue=timepoint_column, palette=palette, ax=ax)
    
    # Add regression lines with confidence intervals for each timepoint
    legend = []
    updated_handles = []  # Store legend handles
    i=0
    for timepoint, subset in data.groupby(timepoint_column):  
        n_samples = len(subset)
        legend.append(f"{timepoint} (n={n_samples})")  # Add sample size to label
        updated_handles.append(mlines.Line2D([], [], color=palette[i], marker='o', linestyle='None', markersize=8))
        i+=1

    sns.regplot(data=data, x=x_column, y=y_column, scatter=False, color='black', ci=95, ax=ax)

    # Labels and title
    ax.set_ylabel(y_column_name, fontsize=fontsize_axes)
    ax.set_title(title, fontsize=fontsize_title)
    ax.grid(True)
    ax.set_xlabel(x_column_name, fontsize=fontsize_axes)
    
    # Adjust legend
    ax.legend_.remove()
    if loc_position != None:
        loc_position = loc_position
        ax.legend(updated_handles, legend, title='', loc=loc_position, fontsize=12, title_fontsize=13)
    
    # Save figure if a path is provided
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    
    return ax


def scatterplot_age(data, x_column, y_column, timepoint_column, title, y_column_name, x_column_name, palette, 
                    ax=None, fontsize_title=16, fontsize_axes=14, save_path=None, loc_position='lower right'):
    """
    Creates a scatter plot with regression lines and confidence intervals for each group in timepoint_column. Also annotates the legend with sample sizes.
    
    Parameters
    ----------
    data : pandas.DataFrame
        The dataset containing values for the x-axis, y-axis, and groups defined by `timepoint_column`.
    x_column : str
        Name of the column in `data` to be plotted on the x-axis.
    y_column : str
        Name of the column in `data` to be plotted on the y-axis.
    timepoint_column : str
        Name of the column that defines groups (e.g., timepoints) used for coloring and legend entries.
    title : str
        Title of the plot.
    y_column_name : str
        Display name (label) for the y-axis.
    x_column_name : str
        Display name (label) for the x-axis.
    palette : list or dict
        Color palette used for different timepoints. Should have at least as many entries as unique groups in `timepoint_column`.
    ax : matplotlib.axes.Axes, optional
        A Matplotlib Axes object to draw the plot on. If None, a new figure and axes are created.
    fontsize_title : int, optional (default=16)
        Font size for the plot title.
    fontsize_axes : int, optional (default=14)
        Font size for the axis labels.
    save_path : str or None, optional
        If provided, the figure will be saved to this file path.
    loc_position : str or None, optional (default='lower right')
        Position of the legend. If None, the legend is omitted.
    """

    data = data.sort_values([timepoint_column])

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter plot
    sns.scatterplot(data=data, x=x_column, y=y_column, hue=timepoint_column, palette=palette, ax=ax)

    # Regression lines per timepoint
    legend = []
    updated_handles = []

    for i, (timepoint, subset) in enumerate(data.groupby(timepoint_column)):
        n_samples = len(subset)
        legend.append(f"{timepoint} (n={n_samples})")
        updated_handles.append(mlines.Line2D([], [], color=palette[i], marker='o', linestyle='None', markersize=8))

        sns.regplot(
            data=subset,
            x=x_column,
            y=y_column,
            scatter=False,
            color=palette[i],
            ci=95,
            ax=ax
        )

    # Labels and title
    ax.set_xlabel(x_column_name, fontsize=fontsize_axes)
    ax.set_ylabel(y_column_name, fontsize=fontsize_axes)
    ax.set_title(title, fontsize=fontsize_title)
    ax.grid(True)

    # Adjust legend
    ax.legend_.remove()
    if loc_position is not None:
        ax.legend(updated_handles, legend, title='', loc=loc_position, fontsize=12, title_fontsize=13)

    # Save figure if a path is provided
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)

    return ax


def plot_estimates(models, titles, x_pos, save_path):
    """
    Creates a multi-panel forest plot for rhythmicity models.

    Parameters:
    -----------
    models : list
        List of fitted model objects.
    titles : list of str
        Titles for each subplot.
    x_pos : list of float
        X-axis positions for the p-value annotations in each subplot.
    figures_path : str
        Directory path where the PDF will be saved.
    """

    predictors = models[0].params.drop(['Intercept', 'Group Var'], errors='ignore').rename(index=name_mapping).index

    fig, axes = plt.subplots(1, len(models), figsize=(22, len(predictors)), sharey=True)

    # Ensure axes is iterable if only one model is provided
    if len(models) == 1:
        axes = [axes]

    for i, (model, ax, title) in enumerate(zip(models, axes, titles)):
        params = model.params.drop(['Intercept', 'Group Var'], errors='ignore').rename(index=name_mapping)
        conf = model.conf_int().drop(['Intercept', 'Group Var'], errors='ignore').rename(index=name_mapping)
        p_values = model.pvalues.drop(['Intercept', 'Group Var'], errors='ignore').rename(index=name_mapping)

        # Align predictors to a consistent order
        params = params.reindex(predictors)
        conf = conf.reindex(predictors)
        p_values = p_values.reindex(predictors)

        # Compute error bars
        lower = params - conf[0]
        upper = conf[1] - params

        # Plot
        ax.errorbar(params.values,
                    params.index,
                    xerr=[lower.values, upper.values],
                    fmt='o',
                    color='black',
                    ecolor='gray',
                    capsize=3)

        ax.axvline(x=0, color='blue', linestyle='--')
        ax.set_title(title, fontsize=16)
        ax.set_xlabel("β Estimate", fontsize=14)
        ax.set_ylabel("")
        ax.set_yticklabels([])
        ax.grid(axis='y')
        ax.set_ylim(-0.5, len(predictors) - 0.5)
        for tick in ax.get_yticklabels():
                    if tick.get_text() in bold_labels:
                        tick.set_fontweight('bold')

        # Add p-values as text annotations next to the points
        pos = x_pos[i]
        for j, (coef, p_val) in enumerate(zip(params.values, p_values.values)):
            if np.isnan(p_val):
                label = "p = nan"
            else:
                label = f"p = {p_val:.3f}" if p_val >= 0.001 else "p < 0.001"
            ax.text(pos, j + 0.1, label, verticalalignment='center', fontsize=10)

    # Label only the first axis
    axes[0].set_yticks(range(len(predictors)))
    axes[0].set_yticklabels(predictors, fontsize=14)
    axes[0].set_ylabel('Predictors', fontsize=14)

    plt.savefig(f"{save_path}", dpi=300, bbox_inches='tight')
    plt.show()


def plot_estimates_vertically(models, titles, save_path):
    """
    Creates a multi-panel forest plot for rhythmicity models.

    Parameters:
    -----------
    models : list
        List of fitted model objects.
    titles : list of str
        Titles for each subplot.
    save_path : str
        Directory path where the PDF will be saved.
    """

    n_models = len(models)
    if n_models <= 2:
        x_size = 5
    else:
        x_size = 12

    fig, axes = plt.subplots(n_models, 1, figsize=(x_size, 5 * n_models), sharex=True)

    # Ensure axes is iterable if only one model is provided
    if n_models == 1:
        axes = [axes]
        x_pos = 50
    else:
        x_pos = 0.3

    for i, (model, ax, title) in enumerate(zip(models, axes, titles)):
        params = model.params.drop(['Intercept', 'Group Var'], errors='ignore').rename(index=name_mapping)
        conf = model.conf_int().drop(['Intercept', 'Group Var'], errors='ignore').rename(index=name_mapping)
        p_values = model.pvalues.drop(['Intercept', 'Group Var'], errors='ignore').rename(index=name_mapping)

        # Compute error bars
        lower = params - conf[0]
        upper = conf[1] - params

        # Plot
        ax.errorbar(params.values,
                    params.index,
                    xerr=[lower.values, upper.values],
                    fmt='o',
                    color='black',
                    ecolor='gray',
                    capsize=3)

        ax.axvline(x=0, color='blue', linestyle='--')
        ax.set_title(title, fontsize=16, pad=10)
        ax.set_xlabel("β Estimate", fontsize=14)
        ax.set_ylabel("Predictors", fontsize=14)
        ax.grid(axis='y')

        # Show predictor names for each model
        ax.set_yticks(range(len(params.index)))
        ax.set_yticklabels(params.index, fontsize=12)
        ax.set_ylim(-0.5, len(params.index) - 0.5)
        
        for tick in ax.get_yticklabels():
            if tick.get_text() in bold_labels:
                tick.set_fontweight('bold')

        # Add p-values as text annotations next to the points
        for j, (coef, p_val) in enumerate(zip(params.values, p_values.values)):
            ax.text(x_pos, j + 0.2, f"p = {p_val:.3f}" if p_val >= 0.001 else "p < 0.001",
                    verticalalignment='center', fontsize=10)

    # Adjust spacing
    plt.subplots_adjust(hspace=0.5)

    # Save and show using the passed save_path
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def prepare_features(feature_table, data_analysis, data_timepoint_sleep, column='babySQUID'):
    """
    Prepare feature and label arrays for model training by merging feature data,
    metadata, and target sleep-quality measurements.
    This function performs the following steps:
    1. Merges `feature_table` with `data_analysis` to add `infant_id` and `timepoint`
       based on their shared index.
    2. Groups features by `infant_id` and `timepoint` and computes the mean for each group.
    3. Merges the aggregated features with sleep-quality data (`data_timepoint_sleep`)
       using the specified target column.
    4. Extracts feature matrix `X`, target vector `y`, and group labels for grouped
       cross-validation.

    Parameters
    ----------
    feature_table : pandas.DataFrame
        A dataframe containing extracted features with an index matching that of 
        `data_analysis`.
    data_analysis : pandas.DataFrame
        Contains at least the columns `infant_id` and `timepoint`, aligned with 
        `feature_table` via index.
    data_timepoint_sleep : pandas.DataFrame
        Table containing sleep-related outcomes, including the specified `column`
        (default: 'babySQUID'), along with `infant_id` and `timepoint`.
    column : str, default='babySQUID'
        The name of the target column to extract from `data_timepoint_sleep`.
    """
    # Merge feature_table with data_analysis on index and specified columns
    feature_table_temp = feature_table.merge(
        data_analysis[['infant_id', 'timepoint']],
        left_index=True,
        right_index=True
    )
    
    # Group by infant_id and timepoint, then compute mean
    feature_table_mean = feature_table_temp.groupby(['infant_id', 'timepoint']).mean()
    
    # Merge with sleep quality data
    feature_table_all = feature_table_mean.merge(
        data_timepoint_sleep[['infant_id', 'timepoint', column]],
        left_index=True,
        right_on=['infant_id', 'timepoint'],
        how='inner'
    )
    
    # Set and reset index
    feature_table_all.set_index(['infant_id', 'timepoint'], inplace=True)
    feature_table_all.reset_index(inplace=True)
    
    # Extract features and labels
    features = feature_table_all.drop(columns=[column, 'infant_id', 'timepoint'])
    X = features.values
    y = feature_table_all[column].values
    groups = feature_table_all['infant_id']
    
    return X, y, groups, features


def evaluate_model(X, y, groups, train_test_split, task_type="regression", n_splits=5, test_size=None):
    """
    Evaluate a machine learning model using grouped cross-validation.

    This function supports both regression and classification tasks and allows for
    multiple grouped splitting strategies (GroupShuffleSplit, GroupKFold, 
    LeaveOneGroupOut). It trains a Random Forest model on each split, computes the
    appropriate evaluation metrics, and returns the best-performing model based on
    R² (regression) or accuracy (classification).

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Feature matrix.
    y : array-like of shape (n_samples,)
        Target vector.
    groups : array-like of shape (n_samples,)
        Group labels used for ensuring train/test splits respect group boundaries.
    train_test_split : {"GroupShuffleSplit", "GroupKFold", "LeaveOneGroupOut"}
        The cross-validation strategy to apply.
    task_type : {"regression", "classification"}, default="regression"
        Determines which model and metrics to use.
    n_splits : int, default=5
        Number of folds or splits (ignored for LeaveOneGroupOut).
    test_size : float, optional
        Test set proportion for GroupShuffleSplit. Ignored for other strategies.
    """
    # Select model based on task
    if task_type == "regression":
        model = RandomForestRegressor(n_estimators=500, random_state=42)
    elif task_type == "classification":
        model = RandomForestClassifier(n_estimators=500, random_state=42)
    else:
        raise ValueError(f"Unsupported task type: {task_type}")

    # Choose splitting strategy
    match train_test_split:
        case 'GroupShuffleSplit':
            random_splits = GroupShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=42)
        case 'GroupKFold':
            random_splits = GroupKFold(n_splits=n_splits)
        case 'LeaveOneGroupOut':
            random_splits = LeaveOneGroupOut()
        case _:
            raise ValueError(f"Unsupported split type: {train_test_split}")

    # Score containers
    r2_scores = []
    mse_scores = []
    mae_scores = []
    acc_scores = []
    f1_scores = []

    best_score = -np.inf
    best_model = None
    best_y_pred = None
    best_y_test = None
    best_X_test = None

    for i, (train_index, test_index) in enumerate(random_splits.split(X, y, groups)):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        if task_type == "regression":
            if len(y_test) < 2:
                r2 = np.nan
            else:
                r2 = r2_score(y_test, y_pred)
                r2_scores.append(r2)

            mse = mean_squared_error(y_test, y_pred)
            mae = mean_absolute_error(y_test, y_pred)
            mse_scores.append(mse)
            mae_scores.append(mae)

            current_score = r2 if not np.isnan(r2) else -np.inf

        elif task_type == "classification":
            acc = accuracy_score(y_test, y_pred)
            f1 = f1_score(y_test, y_pred, average='weighted')
            acc_scores.append(acc)
            f1_scores.append(f1)
            current_score = acc  # Or use f1

        # Save best model
        if current_score > best_score:
            best_score = current_score
            best_model = model
            best_y_pred = y_pred
            best_y_test = y_test
            best_X_test = X_test

    # Summary
    print("\n--- Cross-Validation Summary ---")
    if task_type == "regression":
        avg_r2 = np.nanmean(r2_scores)
        avg_mse = np.mean(mse_scores)
        avg_mae = np.mean(mae_scores)
        print(f"Average R²: {avg_r2:.3f}")
        print(f"Average MSE: {avg_mse:.3f}")
        print(f"Average MAE: {avg_mae:.3f}")
    else:
        avg_acc = np.mean(acc_scores)
        std_acc = np.std(acc_scores)
        avg_f1 = np.mean(f1_scores)
        std_f1 = np.std(f1_scores)
        print(f"Average Accuracy: {avg_acc:.3f} (±{std_acc:.4f})")
        print(f"Average F1 Score (weighted): {avg_f1:.3f} (±{std_f1:.4f})\n")

    return best_model, best_y_pred, best_y_test, best_X_test