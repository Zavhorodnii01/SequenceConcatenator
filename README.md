# Multi-Format Sequence Concatenator with GUI

## Overview

This project provides a user-friendly Graphical User Interface (GUI) for a backend Python class (`SequenceConcatenator`) designed to process multi-gene sequence data from various common bioinformatics file formats (FASTA, Nexus, GenBank). It allows researchers to easily load sequence data for multiple taxa across several gene files, concatenate their sequences into a single alignment, manage missing data, and view/export the resulting alignment, useful statistics, and partition information for downstream phylogenetic analyses.

The GUI simplifies the workflow, making the powerful data processing capabilities accessible through a point-and-click interface without requiring direct interaction with the Python code.

## Features

**Core Backend Features (`SequenceConcatenator.py`):**

*   **Multi-Format Parsing:** Supports reading sequence data from FASTA, Nexus, and GenBank file formats.
*   **Multi-Gene Handling:** Processes data from an arbitrary number of gene files simultaneously.
*   **Multi-Taxa Management:** Automatically identifies all unique taxa present across all input files.
*   **Sequence Concatenation:** Builds a single, long sequence for each taxon by joining their sequences from each gene in order.
*   **Missing Data Handling:** Automatically fills regions corresponding to missing genes for a taxon with gap characters (`-`).
*   **Length Inconsistency Management:** Includes logic to pad shorter sequences with gaps or truncate longer sequences within a gene to ensure alignment column consistency.
*   **Sequence Type Detection:** Attempts to determine if each gene is DNA or Protein based on character composition.
*   **Alignment Statistics:** Calculates and provides statistics about the resulting concatenated alignment, including total length, number of taxa, total missing data, and gene lengths.
*   **Partition Information:** Generates partition boundaries (start-end positions) and estimated data types for each gene.

**GUI Features (`SequenceConcatenatorApp.py`):**

*   **Intuitive Interface:** A simple two-tab layout (`Input` and `Results`) for a clear workflow.
*   **Flexible Input:** Add genes by loading one or multiple files, or by pasting sequence data manually.
*   **Gene Management:** View a list of loaded genes and remove selected ones before processing.
*   **One-Click Processing:** A single "Concatenate & Analyze" button triggers the backend processing.
*   **Comprehensive Results Display:** View the concatenated sequences (formatted as FASTA), detailed statistics, and partition information (formatted for common phylogenetic software like PAUP* and MrBayes) directly within the app.
*   **Results Export:** Easily export the concatenated sequences to a FASTA file, or the partition information and full concatenated alignment to NEXUS files. Includes options for gene-level and codon-position partitions (for DNA genes).
*   **Reset Functionality:** Clear all loaded data and results to start over.
*   **Built with `tkinter`:** Uses standard Python libraries, requiring no external GUI dependencies.

## Requirements

*   Python 3.8 or higher (due to the use of the walrus operator `:=` in the backend class).
*   The standard Python libraries `tkinter`, `ttk`, `filedialog`, `messagebox`, `re`, `os`. (`tkinter` might require separate installation on some Linux/macOS distributions).
*   The `SequenceConcatenator.py` file must be in the same directory as `SequenceConcatenatorApp.py`, as the GUI directly imports and uses the `SequenceConcatenator` class.

## Installation

1.  Clone this repository:

    ```bash
    git clone https://github.com/Zavhorodnii01/SequenceConcatenator.git
    cd SequenceConcatenator
    ```
2.  Ensure both `SequenceConcatenator.py` and `SequenceConcatenatorApp.py` files are present in the same directory.

## Usage

1.  Run the main application script from your terminal:

    ```bash
    python SequenceConcatenatorApp.py
    ```
2.  The GUI window will open.
3.  Navigate to the "Input" tab (this is the default tab).
4.  To add gene data, click the "Add Gene (Manual/File)" button. This will open a new dialog window.
    *   **Manual Input:** In the dialog, you can enter a name for the gene and paste sequence data (preferably in multi-FASTA format) directly into the text area. Click "Add Manual Gene".
    *   **File Loading:** Click "Browse and Load File(s)" in the dialog to open a file chooser. Select one or more gene files (FASTA, Nexus, or GenBank formats are supported). You can optionally provide a base name, or the app will use the filenames. Click "Open" or "Select" and the app will process the file paths.
5.  The names of successfully loaded genes will appear in the "Loaded Genes" list in the main window.
6.  (Optional) To remove genes, select one or more gene names in the list (use Ctrl/Cmd + Click for multiple selections) and click "Remove Selected Gene(s)".
7.  (Optional) Click "Reset All" to clear the list of loaded genes and any previous results.
8.  Once all desired gene files have been added, click the "Concatenate & Analyze" button on the Input tab.
9.  The application will process the data and automatically switch to the "Results" tab.
10. On the "Results" tab, you can view:
    *   The concatenated sequences in FASTA format.
    *   A summary of statistics about the alignment.
    *   Partition information formatted for use in phylogenetic software (examples for PAUP* and MrBayes are provided).
11. Use the "Export FASTA", "Export Partition (NEXUS block)", or "Export Full NEXUS" buttons at the bottom of the Results tab to save the displayed information to files on your computer.

## Output

The application displays and allows export of three main outputs in the "Results" tab:

1.  **Concatenated Sequences:** The combined sequence for each taxon across all processed genes, presented in FASTA format. Missing data is represented by hyphens.
2.  **Statistics:** A summary including the number of taxa, total alignment length, overall percentage of missing data, and missing data counts per taxon and per gene.
3.  **Partition Data:** Information defining the boundaries of each gene within the concatenated alignment. This is provided in a format (NEXUS charset and link blocks) commonly used by phylogenetic software to specify partitions for applying different evolutionary models. Codon-position partitions are also generated for identified DNA genes.
