import os
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kruskal
import scikit_posthocs as sp
import pandas as pd


def process_tetraploid(triads_file, tetraploid_file, diad_file):
    expression_dict = {}
    with open(tetraploid_file, 'r') as tetra_file:
        reader = csv.reader(tetra_file, delimiter='\t')
        next(reader)
        for row in reader:
            key = row[0]
            value = float(row[-1])
            expression_dict[key] = value

    output_data = []
    with open(triads_file, 'r') as triads_file:
        reader = csv.reader(triads_file, delimiter='\t')
        header = next(reader)
        new_header = [
            header[1], f"{header[1]}_expression_values",
            header[2], f"{header[2]}_expression_values",
            "A/G_ratio", "A/G_ratio_log2"
        ]
        output_data.append(new_header)

        for row in reader:
            At_value = expression_dict.get(row[1], "NA")
            G_value = expression_dict.get(row[2], "NA")
            if At_value != "NA" and G_value != "NA" and At_value >= 0.5 and G_value >= 0.5:
                A_G_ratio = At_value / G_value if G_value != 0 else "NaN"
                A_G_ratio_log2 = math.log2(A_G_ratio) if A_G_ratio != "NaN" and A_G_ratio > 0 else "NaN"
            else:
                A_G_ratio = "NaN"
                A_G_ratio_log2 = "NaN"
            output_data.append([row[1], At_value, row[2], G_value, A_G_ratio, A_G_ratio_log2])

    with open(diad_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_data)


def process_files(triads_file, hexaploid_file, triad_file):
    expression_dict = {}
    with open(hexaploid_file, 'r') as hex_file:
        reader = csv.reader(hex_file, delimiter='\t')
        next(reader)
        for row in reader:
            key = row[0]
            value = float(row[-1])
            expression_dict[key] = value

    output_data = []
    with open(triads_file, 'r') as triads_file:
        reader = csv.reader(triads_file, delimiter='\t')
        header = next(reader)
        new_header = [
            header[0], f"{header[0]}_expression_values",
            header[1], f"{header[1]}_expression_values",
            header[2], f"{header[2]}_expression_values",
            "A/G_ratio", "A/G_ratio_log2"
        ]
        output_data.append(new_header)

        for row in reader:
            Am_value = expression_dict.get(row[0], "NA")
            At_value = expression_dict.get(row[1], "NA")
            G_value = expression_dict.get(row[2], "NA")
            if Am_value != "NA" and At_value != "NA" and G_value != "NA":
                A_value = Am_value + At_value
                if A_value >= 0.5 and G_value >= 0.5:
                    A_G_ratio = A_value / G_value if G_value != 0 else "NaN"
                    A_G_ratio_log2 = math.log2(A_G_ratio) if A_G_ratio != "NaN" and A_G_ratio > 0 else "NaN"
                else:
                    A_G_ratio = "NaN"
                    A_G_ratio_log2 = "NaN"
            else:
                A_G_ratio = "NaN"
                A_G_ratio_log2 = "NaN"
            output_data.append([
                row[0], Am_value,
                row[1], At_value,
                row[2], G_value,
                A_G_ratio, A_G_ratio_log2
            ])

    with open(triad_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(output_data)


def merge_datasets(diad_file, triad_file, intersect_file):
    diad_data = []
    with open(diad_file, 'r') as diad:
        reader = csv.DictReader(diad, delimiter='\t')
        for row in reader:
            if row["A/G_ratio"] != "NaN":
                diad_data.append(row)

    triad_data = []
    with open(triad_file, 'r') as triad:
        reader = csv.DictReader(triad, delimiter='\t')
        for row in reader:
            if row["A/G_ratio"] != "NaN":
                triad_data.append(row)

    merged_data = []
    for diad_row in diad_data:
        diad_ag = diad_row["TA2804v2_AG"]
        for triad_row in triad_data:
            if triad_row["TA2804v2_AG"] == diad_ag:
                merged_row = {
                    **diad_row,
                    **{f"triad_{key}": value for key, value in triad_row.items()}
                }
                merged_data.append(merged_row)

    if merged_data:
        with open(intersect_file, 'w', newline='') as output:
            fieldnames = list(diad_data[0].keys()) + [f"triad_{key}" for key in triad_data[0].keys()]
            writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(merged_data)
        print(f"Merged data saved to: {intersect_file}")
    else:
        print("No matching rows found. Merged dataset is empty.")


def calculate_statistics(diad_file, triad_file, stats_output):
    tetraploid = {}
    with open(diad_file, 'r') as diad:
        reader = csv.DictReader(diad, delimiter='\t')
        for row in reader:
            key = f"{row['TA2804v2_AG']}_{row['TA2804v2_GG']}"
            tetraploid[key] = [
                float(row["TA2804v2_AG_expression_values"]),
                float(row["TA2804v2_GG_expression_values"])
            ]

    hexaploid = {}
    with open(triad_file, 'r') as triad:
        reader = csv.DictReader(triad, delimiter='\t')
        for row in reader:
            key = f"{row['T.monococcum.TA10622']}_{row['TA2804v2_AG']}_{row['TA2804v2_GG']}"
            hexaploid[key] = [
                float(row["T.monococcum.TA10622_expression_values"]),
                float(row["TA2804v2_AG_expression_values"]),
                float(row["TA2804v2_GG_expression_values"])
            ]

    tetraploid_values = [value for values in tetraploid.values() for value in values]
    hexaploid_values = [value for values in hexaploid.values() for value in values]

    # Kruskal-Wallis test
    kruskal_result = kruskal(tetraploid_values, hexaploid_values)

    # Prepare data for Dunn's post hoc test
    combined_values = tetraploid_values + hexaploid_values
    groups = ['Tetraploid'] * len(tetraploid_values) + ['Hexaploid'] * len(hexaploid_values)

    # Convert to DataFrame
    df = pd.DataFrame({
        "values": combined_values,
        "groups": groups
    })

    # Dunn's post hoc test
    dunn_result = sp.posthoc_dunn(df, val_col="values", group_col="groups", p_adjust='bonferroni')

    # Write results to stats_output
    with open(stats_output, 'w') as f:
        f.write(f"Kruskal-Wallis Test: Statistic={kruskal_result.statistic}, p-value={kruskal_result.pvalue}\n")
        f.write("\nDunn's Post Hoc Test Results:\n")
        f.write(dunn_result.to_string())
        f.write("\n")

def plot_and_summarize(intersect_file, output_folder):
    diads_no_log = []
    triads_no_log = []
    diads = []
    triads = []
    normalized_triads = []
    normalized_triads_no_log = []

    with open(intersect_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            try:
                diad_no_log = float(row["A/G_ratio"])
                diad_log = float(row["A/G_ratio_log2"])
                triad_no_log = float(row["triad_A/G_ratio"])
                triad_log = float(row["triad_A/G_ratio_log2"])

                diads_no_log.append(diad_no_log)
                diads.append(diad_log)
                triads_no_log.append(triad_no_log)
                triads.append(triad_log)

                if diad_no_log != 0:
                    normalized_value = triad_no_log / diad_no_log
                    normalized_triads_no_log.append(normalized_value)
                    normalized_triads.append(np.log2(normalized_value))
            except ValueError:
                continue

    # Plot the boxplots
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.boxplot([triads, diads], labels=['Triads', 'Diads'])
    plt.axhline(0, color='red', linestyle='--')
    plt.ylabel('A/G Ratio (log2)')
    plt.title('A/G Ratios for Triads and Diads')

    plt.subplot(1, 2, 2)
    plt.boxplot([normalized_triads], labels=['Normalized Triads (log2)'])
    plt.axhline(0, color='red', linestyle='--')
    plt.ylabel('Triad / Diad Ratio (log2)')
    plt.title('Normalized Triads (Log2)')

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, "plots.pdf"))
    plt.close()

    # Calculate summary statistics
    triads_summary = {
        'mean': np.nanmean(triads_no_log),
        'median': np.nanmedian(triads_no_log),
        'std_dev': np.nanstd(triads_no_log)
    }
    diads_summary = {
        'mean': np.nanmean(diads_no_log),
        'median': np.nanmedian(diads_no_log),
        'std_dev': np.nanstd(diads_no_log)
    }
    normalized_triads_no_log_summary = {
        'mean': np.nanmean(normalized_triads_no_log),
        'median': np.nanmedian(normalized_triads_no_log),
        'std_dev': np.nanstd(normalized_triads_no_log)
    }
    normalized_triads_summary = {
        'mean': np.nanmean(normalized_triads),
        'median': np.nanmedian(normalized_triads),
        'std_dev': np.nanstd(normalized_triads)
    }

    # Write summary to a file
    summary_file = os.path.join(output_folder, "summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Summary of A/G Ratios (without log2 conversion):\n\n")
        f.write("Triads:\n")
        f.write(f"  Mean: {triads_summary['mean']:.2f}\n")
        f.write(f"  Median: {triads_summary['median']:.2f}\n")
        f.write(f"  Standard Deviation: {triads_summary['std_dev']:.2f}\n\n")
        f.write("Diads:\n")
        f.write(f"  Mean: {diads_summary['mean']:.2f}\n")
        f.write(f"  Median: {diads_summary['median']:.2f}\n")
        f.write(f"  Standard Deviation: {diads_summary['std_dev']:.2f}\n\n")
        f.write("Normalized Triads (not log2):\n")
        f.write(f"  Mean: {normalized_triads_no_log_summary['mean']:.2f}\n")
        f.write(f"  Median: {normalized_triads_no_log_summary['median']:.2f}\n")
        f.write(f"  Standard Deviation: {normalized_triads_no_log_summary['std_dev']:.2f}\n\n")
        f.write("Normalized Triads (log2):\n")
        f.write(f"  Mean: {normalized_triads_summary['mean']:.2f}\n")
        f.write(f"  Median: {normalized_triads_summary['median']:.2f}\n")
        f.write(f"  Standard Deviation: {normalized_triads_summary['std_dev']:.2f}\n")

    print(f"Summary saved to: {summary_file}")

def main(triads_file, tetraploid_file, hexaploid_file, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    diad_file = os.path.join(output_folder, "diads_ratios.tsv")
    triad_file = os.path.join(output_folder, "triads_ratios.tsv")
    intersect_file = os.path.join(output_folder, "intersect.tsv")
    stats_output = os.path.join(output_folder, "statistics.txt")

    process_tetraploid(triads_file, tetraploid_file, diad_file)
    process_files(triads_file, hexaploid_file, triad_file)
    merge_datasets(diad_file, triad_file, intersect_file)
    calculate_statistics(diad_file, triad_file, stats_output)
    plot_and_summarize(intersect_file, output_folder)


# Run the pipeline
triads_file = "/path/to/triads/file"
tetraploid_file = "/path/to/expression/file"
hexaploid_file = "/path/to/expression/file"
output_folder = f"gene_dosage_{os.path.basename(tetraploid_file).split('.')[0]}_{os.path.basename(hexaploid_file).split('.')[0]}"
main(triads_file, tetraploid_file, hexaploid_file, output_folder)

