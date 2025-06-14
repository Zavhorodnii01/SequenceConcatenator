import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import re
import os

from SequenceConcatenator import SequenceConcatenator
class SequenceConcatenatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Sequence Concatenator UI")
        self.root.geometry("1000x800") # Increased size for more content
        self.root.configure(bg="#f4f4f4")


        self.gene_names = []
        self.raw_gene_contents = [] # Stores raw content lines for each gene file/input

        # --- Styling ---
        self.bg_color = "#f4f4f4"
        self.primary_color = "#1e3a8a" # blue
        self.secondary_color = "#888888" # grey
        self.text_color = "black"
        self.heading_color = "#1e3a8a" # blue
        self.button_fg = "white"
        self.button_bg_primary = self.primary_color
        self.button_bg_secondary = self.secondary_color
        self.active_bg_primary = "#144c8f"
        self.active_bg_secondary = "#555555"
        self.font_normal = ("Arial", 11)
        self.font_bold = ("Arial", 11, "bold")
        self.font_heading = ("Arial", 12, "bold")
        self.font_button = ("Segoe UI", 10, "bold")
        self.font_mono = ("Courier New", 10) # Monospaced for sequence display

        # Notebook (Tabs)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(pady=10, padx=10, fill="both", expand=True)

        # Create tabs
        # Using ttk.Frame for consistent styling with the notebook
        self.tabs = {
            "Input": ttk.Frame(self.notebook, padding="10", style='TFrame'),
            "Results": ttk.Frame(self.notebook, padding="10", style='TFrame'), # New tab for results
        }

        # Configure ttk Style
        style = ttk.Style()
        style.theme_use('default') # Use default theme, customize from there
        style.configure('TFrame', background=self.bg_color)
        style.configure('TNotebook.Tab', padding=[10, 5], font=self.font_normal)
        style.map('TNotebook.Tab',
                  background=[('selected', self.primary_color)],
                  foreground=[('selected', 'white')],
                  expand=[('selected', [1, 1, 1, 0])]) # Keep selected tab slightly expanded

        for tab_name, tab_frame in self.tabs.items():
            self.notebook.add(tab_frame, text=tab_name)

        self.setup_input_tab()


        self.results_content_frame = ttk.Frame(self.tabs["Results"], style='TFrame')
        self.results_content_frame.pack(fill="both", expand=True)
        self.setup_results_tab(self.results_content_frame)


        # Store results internally for export
        self._concatenated_sequences = None
        # Now expects partition data in the format: [(gene_name, 'start-end', gene_type), ...]
        self._partition_data = None
        self._statistics = None


    def setup_input_tab(self):
        input_frame = self.tabs["Input"]
        # Configure grid layout for better control
        input_frame.columnconfigure(0, weight=1) # Allow the list container to expand
        input_frame.rowconfigure(0, weight=1) # Allow the list container to expand

        # Area to display loaded gene names
        list_container_frame = tk.Frame(input_frame, bg=self.bg_color)
        list_container_frame.grid(row=0, column=0, pady=10, padx=10, sticky="nsew") # Use grid

        tk.Label(list_container_frame, text="Loaded Genes:", font=self.font_heading, bg=self.bg_color, fg=self.heading_color).pack(anchor="w", pady=(0, 5))

        # Frame for Listbox and Scrollbar
        listbox_frame = tk.Frame(list_container_frame, bg="white", bd=1, relief="sunken")
        listbox_frame.pack(fill="both", expand=True)

        self.gene_listbox = tk.Listbox(
            listbox_frame,
            font=self.font_normal,
            fg=self.text_color,
            bg="white",
            selectbackground=self.primary_color,
            selectforeground="white",
            activestyle="none",
            bd=0,
            highlightthickness=0,
            exportselection=False,
            selectmode=tk.EXTENDED # Allow multiple selections
        )
        self.gene_listbox.pack(side="left", fill="both", expand=True)

        list_scrollbar = tk.Scrollbar(listbox_frame, orient="vertical", command=self.gene_listbox.yview)
        list_scrollbar.pack(side="right", fill="y")
        self.gene_listbox.config(yscrollcommand=list_scrollbar.set)


        # Buttons at the bottom
        btn_frame = tk.Frame(input_frame, bg=self.bg_color)
        btn_frame.grid(row=1, column=0, sticky="ew", pady=10, padx=10) # Use grid

        # Use grid for buttons within btn_frame for better alignment
        btn_frame.columnconfigure(0, weight=1) # Allow left buttons to push right buttons
        btn_frame.columnconfigure(1, weight=1)
        btn_frame.columnconfigure(2, weight=0) # Don't expand action buttons
        btn_frame.columnconfigure(3, weight=0)


        # "Add Gene (Manual/File)" button - This button opens the dialog
        add_gene_btn = tk.Button(
            btn_frame,
            text="Add Gene (Manual/File)",
            command=self.open_add_gene_dialog,
            bg=self.button_bg_primary,
            fg=self.button_fg,
            activebackground=self.active_bg_primary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        add_gene_btn.grid(row=0, column=0, sticky="w", padx=(0, 10)) # Use grid

        # "Remove Selected Gene" button
        remove_gene_btn = tk.Button(
            btn_frame,
            text="Remove Selected Gene(s)",
            command=self.remove_selected_gene,
            bg=self.button_bg_secondary,
            fg=self.button_fg,
            activebackground=self.active_bg_secondary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        remove_gene_btn.grid(row=0, column=1, sticky="w", padx=(0, 10)) # Use grid


        # Right-aligned Submit and Reset buttons - put them in a sub-frame for alignment
        action_btn_frame = tk.Frame(btn_frame, bg=self.bg_color)
        action_btn_frame.grid(row=0, column=3, sticky="e") # Put action frame on the right

        reset_btn = tk.Button(
            action_btn_frame,
            text="Reset All",
            command=self.on_reset,
            bg=self.button_bg_secondary,
            fg=self.button_fg,
            activebackground=self.active_bg_secondary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        reset_btn.pack(side="left", padx=(0, 10))

        submit_btn = tk.Button(
            action_btn_frame,
            text="Concatenate & Analyze",
            command=self.on_submit,
            bg=self.button_bg_primary,
            fg=self.button_fg,
            activebackground=self.active_bg_primary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        submit_btn.pack(side="left")


    def setup_results_tab(self, parent_frame):
        # results_frame = self.tabs["Results"] # We now use the parent_frame passed in
        parent_frame.columnconfigure(0, weight=1)
        parent_frame.rowconfigure(0, weight=1) # Statistics/Partition row
        parent_frame.rowconfigure(1, weight=2) # Sequences row
        parent_frame.rowconfigure(2, weight=0) # Buttons row


        # Frame for Statistics and Partition
        # Use a single frame with packed sub-frames inside, placed in grid row 0
        top_panes_frame = ttk.Frame(parent_frame, style='TFrame')
        top_panes_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        top_panes_frame.columnconfigure(0, weight=1)
        top_panes_frame.columnconfigure(1, weight=1)
        top_panes_frame.rowconfigure(0, weight=1)


        # Frame for Statistics Display ---
        stats_frame = tk.LabelFrame(top_panes_frame, text="Statistics", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=10)
        stats_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5) # Use grid
        stats_frame.columnconfigure(0, weight=1)
        stats_frame.rowconfigure(0, weight=1)


        self.stats_text = tk.Text(
            stats_frame,
            font=self.font_normal,
            bg="white",
            fg=self.text_color,
            height=8, # Give it a starting height
            state="disabled",
            wrap="word" # word wrap is fine for stats
        )
        self.stats_text.grid(row=0, column=0, sticky="nsew")


        # Frame for Partition Display ---
        partition_frame = tk.LabelFrame(top_panes_frame, text="Partition (NEXUS format example)", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=10)
        partition_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5) # Use grid
        partition_frame.columnconfigure(0, weight=1)
        partition_frame.rowconfigure(0, weight=1)


        partition_text_frame = tk.Frame(partition_frame) # Container for text and scrollbars
        partition_text_frame.grid(row=0, column=0, sticky="nsew")
        partition_text_frame.columnconfigure(0, weight=1)
        partition_text_frame.rowconfigure(0, weight=1)


        self.partition_text = tk.Text(
            partition_text_frame,
            font=self.font_mono,
            bg="white",
            fg=self.text_color,
            height=8, # Give it a starting height
            state="disabled",
            wrap="none" # Use horizontal scrollbar for partition display
        )
        self.partition_text.grid(row=0, column=0, sticky="nsew")

        # Add scrollbars to the partition text frame
        partition_scrollbar_y = tk.Scrollbar(partition_text_frame, orient="vertical", command=self.partition_text.yview)
        partition_scrollbar_y.grid(row=0, column=1, sticky="ns")
        self.partition_text.config(yscrollcommand=partition_scrollbar_y.set)

        partition_scrollbar_x = tk.Scrollbar(partition_frame, orient="horizontal", command=self.partition_text.xview) # Scrollbar at the bottom of the labelframe
        partition_scrollbar_x.grid(row=1, column=0, sticky="ew")
        self.partition_text.config(xscrollcommand=partition_scrollbar_x.set)


        # --- Concatenated Sequences Display ---
        # Placed in grid row 1
        sequences_frame = tk.LabelFrame(parent_frame, text="Concatenated Sequences (FASTA Format)", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=10)
        sequences_frame.grid(row=1, column=0, sticky="nsew", padx=5, pady=5) # Use grid
        sequences_frame.columnconfigure(0, weight=1)
        sequences_frame.rowconfigure(0, weight=1)

        seq_text_frame = tk.Frame(sequences_frame) # Container for text and scrollbars
        seq_text_frame.grid(row=0, column=0, sticky="nsew")
        seq_text_frame.columnconfigure(0, weight=1)
        seq_text_frame.rowconfigure(0, weight=1)

        self.concatenated_seqs_text = tk.Text(
            seq_text_frame,
            font=self.font_mono,
            bg="white",
            fg=self.text_color,
            state="disabled",
            wrap="none" # Use horizontal scrollbar for sequences
        )
        self.concatenated_seqs_text.grid(row=0, column=0, sticky="nsew")

        # Add scrollbars to the sequence text frame
        seq_scrollbar_y = tk.Scrollbar(seq_text_frame, orient="vertical", command=self.concatenated_seqs_text.yview)
        seq_scrollbar_y.grid(row=0, column=1, sticky="ns")
        self.concatenated_seqs_text.config(yscrollcommand=seq_scrollbar_y.set)

        seq_scrollbar_x = tk.Scrollbar(sequences_frame, orient="horizontal", command=self.concatenated_seqs_text.xview) # Scrollbar at the bottom of the labelframe
        seq_scrollbar_x.grid(row=1, column=0, sticky="ew")
        self.concatenated_seqs_text.config(xscrollcommand=seq_scrollbar_x.set)


        # --- Export Buttons ---
        # Placed in grid row 2
        export_btn_frame = tk.Frame(parent_frame, bg=self.bg_color)
        export_btn_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10), padx=10) # Use grid

        # Use pack within the export frame for simple arrangement
        export_fasta_btn = tk.Button(
            export_btn_frame,
            text="Export FASTA",
            command=self.export_fasta,
            bg=self.button_bg_primary,
            fg=self.button_fg,
            activebackground=self.active_bg_primary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_fasta_btn.pack(side="left", padx=(0, 10))

        export_partition_btn = tk.Button(
            export_btn_frame,
            text="Export Partition (NEXUS block)",
            command=self.export_partition,
            bg=self.button_bg_primary,
            fg=self.button_fg,
            activebackground=self.active_bg_primary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_partition_btn.pack(side="left", padx=(0, 10))

        # Added button for full NEXUS export including concatenated sequences
        export_full_nexus_btn = tk.Button(
             export_btn_frame,
             text="Export Full NEXUS",
             command=self.export_full_nexus,
             bg=self.button_bg_primary,
             fg=self.button_fg,
             activebackground=self.active_bg_primary,
             activeforeground=self.button_fg,
             font=self.font_button,
             bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_full_nexus_btn.pack(side="left", padx=(0, 10))


    def open_add_gene_dialog(self):
        """Opens a new window to add a gene manually or load from file."""
        dialog = tk.Toplevel(self.root)
        dialog.title("Add New Gene or Load File")
        dialog.transient(self.root) # Keep dialog on top of main window
        dialog.grab_set() # Modal dialog - prevents interaction with main window
        dialog.configure(bg=self.bg_color)
        # Calculate a reasonable size and center the dialog relative to the main window
        main_x = self.root.winfo_x()
        main_y = self.root.winfo_y()
        main_width = self.root.winfo_width()
        main_height = self.root.winfo_height()
        dialog_width = 600 # Increased width slightly
        dialog_height = 700 # Increased height significantly to make space for text area
        dialog_x = main_x + (main_width // 2) - (dialog_width // 2)
        dialog_y = main_y + (main_height // 2) - (dialog_height // 2)
        dialog.geometry(f"{dialog_width}x{dialog_height}+{dialog_x}+{dialog_y}")
        dialog.resizable(True, True) # Allow resizing


        # Use grid for the dialog frame
        dialog_frame = tk.Frame(dialog, bg=self.bg_color)
        dialog_frame.pack(fill="both", expand=True, padx=15, pady=10) # Use pack for the main dialog_frame

        # Configure grid for contents *within* dialog_frame
        dialog_frame.columnconfigure(0, weight=1)
        # Manual section (row 0) should take most vertical space
        dialog_frame.rowconfigure(0, weight=3) # Give manual section more weight
        # File section (row 1) fixed size (weight=0)
        dialog_frame.rowconfigure(1, weight=0)
        # Button section (row 2) fixed size (weight=0)
        dialog_frame.rowconfigure(2, weight=0)


        # --- Manual Input Section ---
        manual_frame = tk.LabelFrame(dialog_frame, text="Add Manually (Paste Multi-FASTA Data)", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=10)
        manual_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 10)) # Use grid
        # Configure grid for contents *within* manual_frame
        manual_frame.columnconfigure(0, weight=1)
        manual_frame.rowconfigure(0, weight=0) # Label - fixed height
        manual_frame.rowconfigure(1, weight=0) # Entry - fixed height
        manual_frame.rowconfigure(2, weight=0) # Label - fixed height
        manual_frame.rowconfigure(3, weight=1) # Text frame - EXPANDS
        manual_frame.rowconfigure(4, weight=0) # Button - fixed height


        tk.Label(manual_frame, text="Gene Name :", font=self.font_normal, bg=self.bg_color, fg=self.text_color).grid(row=0, column=0, sticky="w")
        name_entry = tk.Entry(manual_frame, font=self.font_normal)
        name_entry.grid(row=1, column=0, sticky="ew", pady=(0, 5))

        tk.Label(manual_frame, text="Sequence Data (Paste Multi-FASTA here):", font=self.font_normal, bg=self.bg_color, fg=self.text_color).grid(row=2, column=0, sticky="nw")

        # Frame for text widget and scrollbar
        manual_text_frame = tk.Frame(manual_frame)
        manual_text_frame.grid(row=3, column=0, sticky="nsew", pady=(0, 5)) # This frame expands
        manual_text_frame.columnconfigure(0, weight=1) # Text widget expands horizontally
        manual_text_frame.rowconfigure(0, weight=1) # Text widget expands vertically


        seq_text = tk.Text(manual_text_frame, font=self.font_mono, wrap="word") # Let manual entry wrap
        seq_text.grid(row=0, column=0, sticky="nsew") # Text widget fills its frame
        manual_scrollbar_y = tk.Scrollbar(manual_text_frame, orient="vertical", command=seq_text.yview)
        manual_scrollbar_y.grid(row=0, column=1, sticky="ns") # Scrollbar fills vertically
        seq_text.config(yscrollcommand=manual_scrollbar_y.set)


        add_manual_btn = tk.Button(
            manual_frame,
            text="Add Manual Gene",
            command=lambda: self.add_gene_from_dialog(dialog, name_entry.get(), seq_text.get("1.0", "end").strip()),
            bg=self.button_bg_primary,
            fg=self.button_fg,
            activebackground=self.active_bg_primary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=10, pady=5, cursor="hand2"
        )
        add_manual_btn.grid(row=4, column=0, pady=(5,0))


        # --- Load File Section ---
        load_frame = tk.LabelFrame(dialog_frame, text="Load from File(s) (FASTA, NEXUS, GenBank)", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=10)
        load_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10)) # Use grid (does not expand vertically)
        # Configure grid for contents *within* load_frame
        load_frame.columnconfigure(0, weight=1)
        load_frame.rowconfigure(0, weight=0) # Label - fixed height
        load_frame.rowconfigure(1, weight=0) # Entry - fixed height
        load_frame.rowconfigure(2, weight=0) # Button - fixed height


        tk.Label(load_frame, text="Gene Name for File(s) :", font=self.font_normal, bg=self.bg_color, fg=self.text_color).grid(row=0, column=0, sticky="w")
        load_name_entry = tk.Entry(load_frame, font=self.font_normal)
        load_name_entry.grid(row=1, column=0, sticky="ew", pady=(0, 5))

        load_file_btn = tk.Button(
            load_frame,
            text="Browse and Load File(s)",
            command=lambda: self.load_file_dialog_and_store_raw(dialog, load_name_entry.get().strip()), # Store raw here
            bg=self.button_bg_primary,
            fg=self.button_fg,
            activebackground=self.active_bg_primary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=10, pady=5, cursor="hand2"
        )
        load_file_btn.grid(row=2, column=0, pady=(5,0))


        # --- Close Button ---
        close_btn = tk.Button(
            dialog_frame,
            text="Cancel",
            command=dialog.destroy,
            bg=self.button_bg_secondary,
            fg=self.button_fg,
            activebackground=self.active_bg_secondary,
            activeforeground=self.button_fg,
            font=self.font_button,
            bd=0, relief="flat", padx=10, pady=5, cursor="hand2"
        )
        close_btn.grid(row=2, column=0, pady=10)


        # Set focus to the dialog
        dialog.focus_set()

    def add_gene_from_dialog(self, dialog, gene_name, raw_content):
        """Adds a single gene from the manual input fields (raw text)."""
        gene_name = gene_name.strip()
        raw_content = raw_content.strip()

        if not gene_name:
            messagebox.showwarning("Warning", "Please provide a gene name.", parent=dialog)
            return
        if not raw_content:
            messagebox.showwarning("Warning", "Please paste sequence data.", parent=dialog)
            return

        # Store the raw content as lines
        raw_lines = raw_content.splitlines()

        # Add the raw gene data
        self._add_raw_gene_data(gene_name, raw_lines) # Pass raw lines
        dialog.destroy() # Close the dialog after adding


    def load_file_dialog_and_store_raw(self, dialog, manual_gene_name=""):
        """Opens file dialog, loads selected file(s), stores raw content."""
        filepaths = filedialog.askopenfilenames(
            initialdir=".", # Could make this user-configurable later
            title="Select Gene File(s)",
            filetypes=(
                ("Sequence files", "*.fasta *.fas *.fna *.faa *.fa *.nexus *.nex *.gbff"),
                ("FASTA files", "*.fasta *.fas *.fna *.faa *.fa"),
                ("NEXUS files", "*.nexus *.nex"),
                ("GenBank files", "*.gbff"),
                ("All files", "*.*")
            ),
            parent=dialog # Make dialog parent
        )
        if not filepaths:
            return # User cancelled

        successfully_loaded_names = []
        failed_files = []

        for filepath in filepaths:
            try:
                with open(filepath, "r", encoding='utf-8', errors='ignore') as f: # Use utf-8, ignore errors for robustness
                    raw_content_lines = f.readlines()

                if not raw_content_lines or all(line.strip() == '' for line in raw_content_lines):
                    messagebox.showwarning("Warning", f"File '{os.path.basename(filepath)}' is empty or contains only whitespace. Skipping.", parent=dialog)
                    failed_files.append(f"{os.path.basename(filepath)} (Empty)")
                    continue

                # Use manual name if provided, otherwise use filename without extension
                # Note: If loading multiple files with a single manual name, all will get the same base name,
                # which _add_raw_gene_data handles by adding _1, _2 etc.
                gene_name = manual_gene_name if manual_gene_name else os.path.splitext(os.path.basename(filepath))[0]

                # Add the raw gene data
                added_name = self._add_raw_gene_data(gene_name, raw_content_lines) # Store raw lines
                successfully_loaded_names.append(added_name)

            except Exception as e:
                messagebox.showerror("File Reading Error", f"Error reading file '{os.path.basename(filepath)}':\n{e}", parent=dialog)
                failed_files.append(f"{os.path.basename(filepath)} (Error: {e})")
                # Continue to try other files

        # Report overall status
        status_msg = ""
        if successfully_loaded_names:
             status_msg += "Successfully loaded the following gene file(s) (parsing happens on submit):\n" + "\n".join(successfully_loaded_names)
        if failed_files:
             if status_msg: status_msg += "\n\n"
             status_msg += "Failed to load the following file(s):\n" + "\n".join(failed_files)

        if status_msg:
             if successfully_loaded_names:
                 # If any files were successfully loaded, show info and close the dialog
                 messagebox.showinfo("Load Status", status_msg, parent=dialog)
                 dialog.destroy()
             else:
                  # If only failures, show error message and leave dialog open for user to try again or manually add
                  messagebox.showerror("Load Failed", status_msg, parent=dialog)
        else:
             # Should not happen if files were selected, and if no files were selected, we returned earlier.
             pass


    def _add_raw_gene_data(self, gene_name, raw_lines):
        """Internal method to add raw gene content (list of lines) and update the listbox."""
        if not gene_name or not raw_lines:

            return

        # Handle potential duplicate names in the UI listbox display
        original_name = gene_name
        counter = 1
        # Check against names already in the listbox display
        existing_names_in_listbox = list(self.gene_listbox.get(0, tk.END))
        temp_name = gene_name
        while temp_name in existing_names_in_listbox:
            temp_name = f"{original_name}_{counter}"
            counter += 1
        gene_name = temp_name # Use the unique name for display and internal storage


        # Store name and raw content lines, keeping lists synchronized
        self.gene_names.append(gene_name)
        self.raw_gene_contents.append(raw_lines) # Store the raw lines
        self.gene_listbox.insert(tk.END, gene_name)

        return gene_name # Return the name actually added (might have _counter)


    def remove_selected_gene(self):
        """Removes the selected gene(s) from the listbox and internal raw data."""
        selected_indices = self.gene_listbox.curselection()
        if not selected_indices:
            messagebox.showwarning("Warning", "Please select one or more genes to remove.")
            return

        # Get the names of the genes being removed for the confirmation message
        selected_gene_names = [self.gene_listbox.get(i) for i in selected_indices]
        confirm = messagebox.askyesno("Confirm Removal", f"Are you sure you want to remove the selected gene(s)?\n\n" + "\n".join(selected_gene_names))
        if not confirm:
            return

        # Remove from the listbox and internal data structures
        # Iterate backwards through selected indices to avoid issues with index changes
        for i in sorted(selected_indices, reverse=True):

            self.gene_listbox.delete(i)
            # Remove from internal lists at the corresponding index
            del self.gene_names[i]
            del self.raw_gene_contents[i] # Remove the raw lines for this gene



        # Also clear results if genes are removed, as results are no longer valid
        self.clear_results_display()


    def on_submit(self):
        """Handles the submit action - passes raw data to backend,
           calls the SequenceConcatenator, and displays results."""
        if not self.raw_gene_contents:
            messagebox.showwarning("Warning", "Please add at least one gene before submitting.")
            return

        # --- Pass raw data to the backend ---
        # The backend SequenceConcatenator's __init__ now directly accepts the list of raw content lists.
        backend_input = self.raw_gene_contents

        self.clear_results_display() # Clear previous results

        try:
            # --- Instantiate and Run SequenceConcatenator ---
            concatenator = SequenceConcatenator(backend_input)

            # --- Get Results ---
            # The backend's get_partition now returns (gene_name, 'start-end', gene_type)
            self._concatenated_sequences = concatenator.get_concatenated_sequences()
            self._statistics = concatenator.get_statistics()
            self._partition_data = concatenator.get_partition()

            if not self._concatenated_sequences:
                 messagebox.showwarning("Concatenation Warning", "No concatenated sequences were produced by the backend. This might happen if no valid data was parsed or there are no common taxa across genes.")
                 self.clear_results_display() # Ensure results tab is empty if concatenation failed
                 return

            # --- Display Results ---
            self.display_results(self._concatenated_sequences, self._statistics, self._partition_data)

            # --- Switch to Results Tab ---
            self.notebook.select(self.tabs["Results"])


        except Exception as e:
            # Handle errors that might occur in your SequenceConcatenator class during parsing/concatenation
            messagebox.showerror("Processing Error", f"An error occurred during sequence processing in the backend.\nError details: {e}")

            import traceback
            traceback.print_exc()
            # Clear any partial results and stored data
            self._concatenated_sequences = None
            self._partition_data = None
            self._statistics = None
            self.clear_results_display()


    def display_results(self, concatenated_sequences, statistics, partition):
        """Populates the Results tab widgets with data."""
        # Enable editing to insert text
        self.concatenated_seqs_text.config(state="normal")
        self.stats_text.config(state="normal")
        self.partition_text.config(state="normal")

        # Clear previous content
        self.concatenated_seqs_text.delete("1.0", tk.END)
        self.stats_text.delete("1.0", tk.END)
        self.partition_text.delete("1.0", tk.END)

        # Display Concatenated Sequences (in FASTA format)
        self.concatenated_seqs_text.insert(tk.END, "Concatenated Sequences (FASTA):\n\n")
        if concatenated_sequences:
            # Sort by taxon name for consistent display
            sorted_taxons = sorted(concatenated_sequences.keys())
            for taxon in sorted_taxons:
                sequence = concatenated_sequences[taxon]

                wrapped_seq = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
                self.concatenated_seqs_text.insert(tk.END, f">{taxon}\n{wrapped_seq}\n")
        else:
            self.concatenated_seqs_text.insert(tk.END, "No concatenated sequences produced by backend.\n")

        # Display Statistics
        self.stats_text.insert(tk.END, "Statistics:\n\n")
        if statistics:
            # Format statistics nicely. Ensure keys are displayed in a readable order.
            # Define a preferred display order for common stats
            preferred_order = [
                "Number of Taxa",
                "Number of Genes",
                "Total Length",
                "Percentage Overall Missing Data (%)",
                "Missing Data Per Taxon",
                "Missing Data Per Gene",
                "Taxon-Gene Matrix Sparsity (%)",
            ]
            displayed_keys = set()

            # Display preferred keys first
            for key in preferred_order:
                 if key in statistics:
                      if isinstance(statistics[key], dict):
                          self.stats_text.insert(tk.END, f"{key}:\n")
                          # Sort dict keys for consistent display
                          for sub_key in sorted(statistics[key].keys()):
                               sub_value = statistics[key][sub_key]
                               # Format float values
                               if isinstance(sub_value, float):
                                    self.stats_text.insert(tk.END, f"  {sub_key}: {sub_value:.2f}\n")
                               else:
                                     self.stats_text.insert(tk.END, f"  {sub_key}: {sub_value}\n")
                      # Format float values at top level
                      elif isinstance(statistics[key], float):
                           self.stats_text.insert(tk.END, f"{key}: {statistics[key]:.2f}\n")
                      else:
                           self.stats_text.insert(tk.END, f"{key}: {statistics[key]}\n")
                      displayed_keys.add(key) # Mark as displayed

            # Display any other stats returned by the backend
            for key, value in statistics.items():
                if key not in displayed_keys: # Avoid displaying duplicates
                     if isinstance(value, dict):
                          self.stats_text.insert(tk.END, f"{key}:\n")
                          for sub_key in sorted(value.keys()): # Sort sub-keys
                               sub_value = value[sub_key]
                               # Format float values
                               if isinstance(sub_value, float):
                                    self.stats_text.insert(tk.END, f"  {sub_key}: {sub_value:.2f}\n")
                               else:
                                     self.stats_text.insert(tk.END, f"  {sub_key}: {sub_value}\n")
                     # Format float values at top level
                     elif isinstance(value, float):
                          self.stats_text.insert(tk.END, f"{key}: {value:.2f}\n")
                     else:
                          self.stats_text.insert(tk.END, f"{key}: {value}\n")
        else:
            self.stats_text.insert(tk.END, "No statistics available from backend.\n")

        # Display Partition (NEXUS format example)
        self.partition_text.insert(tk.END, "Partition (NEXUS format example):\n\n")
        if partition:
             # partition is expected to be [(gene_name, 'start-end', gene_type), ...]
             self.partition_text.insert(tk.END, "#nexus\n\n")
             self.partition_text.insert(tk.END, "BEGIN PAUP;\n")

             self.partition_text.insert(tk.END, "  [ Charsets for gene partitions ]\n")
             for backend_gene_name, coord_range, gene_type in partition:
                 clean_gene_name = re.sub(r'\W+', '_', backend_gene_name) # Clean name for Nexus charset identifier
                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name} = {coord_range};\n")

             # Add Codon Partitions for DNA genes
             self.partition_text.insert(tk.END, "\n  [ Charsets for codon positions (DNA genes only) ]\n")
             dna_genes_with_partitions = []
             for backend_gene_name, coord_range, gene_type in partition:
                  if gene_type.upper() == 'DNA': # Case-insensitive check
                       clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                       try:
                            start, end = map(int, coord_range.split('-'))
                            # Check if the range is long enough for at least one codon (3 bp)
                            # Also handle 1-based vs 0-based. NEXUS is typically 1-based.
                            # The backend range is 1-based: '1-100' means positions 1 through 100 inclusive.
                            # Length is end - start + 1.
                            length = end - start + 1
                            if length >= 3:
                                 # Codon positions are relative to the start of the gene segment in the concatenated sequence
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n")
                                 dna_genes_with_partitions.append(clean_gene_name)
                            elif length > 0:
                                 # Gene too short for full codons, but might have partial data
                                 self.partition_text.insert(tk.END, f"  [ Warning: Gene {clean_gene_name} ({coord_range}) is too short ({length} bp) for full codon partitions. ]\n")
                       except ValueError:
                            self.partition_text.insert(tk.END, f"  [ Warning: Could not parse range '{coord_range}' for codon partitions for gene {clean_gene_name}. ]\n")
                  elif gene_type.upper() == 'PROTEIN':
                       self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")


             if not dna_genes_with_partitions:
                  self.partition_text.insert(tk.END, "  [ No DNA genes found or genes too short for codon partitions. ]\n")


             # Add LINK block
             self.partition_text.insert(tk.END, "\n  [ Link block (gene level) ]\n")
             # Only add LINK block if there's something to link
             if partition:
                 self.partition_text.insert(tk.END, "  link characters = ")
                 link_strings = []
                 for backend_gene_name, coord_range, gene_type in partition:
                      clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                      link_strings.append(f"{clean_gene_name} : {coord_range}")
                 self.partition_text.insert(tk.END, ", ".join(link_strings) + ";\n")
             else:
                  self.partition_text.insert(tk.END, "  [ No genes to link ]\n")



             if dna_genes_with_partitions:
                 self.partition_text.insert(tk.END, "\n  [ Link block (by codon position, DNA genes only) ]\n")
                 self.partition_text.insert(tk.END, "  link characters = ")
                 codon_link_strings = []
                 # Collect pos1, pos2, pos3 ranges from *all* DNA genes that were long enough
                 # Need to recreate the ranges as PAUP link format is `charset_name : range`
                 pos1_ranges = []
                 pos2_ranges = []
                 pos3_ranges = []
                 for name, range_str, gene_type in partition:
                      if gene_type.upper() == 'DNA':
                           clean_name = re.sub(r'\W+', '_', name)
                           try:
                                start, end = map(int, range_str.split('-'))
                                length = end - start + 1
                                if length >= 3:
                                     pos1_ranges.append(f"{clean_name}_pos1 : {start}-{end}\\3")
                                     pos2_ranges.append(f"{clean_name}_pos2 : {start+1}-{end}\\3")
                                     pos3_ranges.append(f"{clean_name}_pos3 : {start+2}-{end}\\3")
                           except ValueError:
                                pass # Ignore errors already warned about

                 all_codon_ranges = pos1_ranges + pos2_ranges + pos3_ranges
                 if all_codon_ranges:
                     self.partition_text.insert(tk.END, ", ".join(all_codon_ranges) + ";\n")
                 else:
                     self.partition_text.insert(tk.END, "[ No DNA genes long enough for codon link block ]\n")


             self.partition_text.insert(tk.END, "END; [PAUP]\n\n")

             # Add MrBayes block example
             self.partition_text.insert(tk.END, "BEGIN mrbayes;\n")
             self.partition_text.insert(tk.END, "  [ Example MrBayes partition block ]\n")
             # Add charsets again for MrBayes block (both gene and codon levels)
             self.partition_text.insert(tk.END, "  [ Gene level charsets ]\n")
             gene_charset_names = [] # For 'partition by_gene' command
             for backend_gene_name, coord_range, gene_type in self._partition_data:
                 clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name} = {coord_range};\n")
                 gene_charset_names.append(clean_gene_name)

             self.partition_text.insert(tk.END, "\n  [ Codon position charsets (DNA genes only) ]\n")
             mr_bayes_dna_charsets_exist = False
             pos1_names = [] # For 'partition by_codon_pos' command
             pos2_names = []
             pos3_names = []

             for backend_gene_name, coord_range, gene_type in self._partition_data:
                  if gene_type.upper() == 'DNA':
                       clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                       try:
                            start, end = map(int, coord_range.split('-'))
                            length = end - start + 1
                            if length >= 3:
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n")
                                 mr_bayes_dna_charsets_exist = True
                                 pos1_names.append(f"{clean_gene_name}_pos1")
                                 pos2_names.append(f"{clean_gene_name}_pos2")
                                 pos3_names.append(f"{clean_gene_name}_pos3")
                       except ValueError:
                               pass # Error message handled in PAUP block display
                  elif gene_type.upper() == 'PROTEIN':
                        self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")


             if not mr_bayes_dna_charsets_exist:
                  self.partition_text.insert(tk.END, "  [ No DNA genes found or genes too short for codon partitions. ]\n")

             self.partition_text.insert(tk.END, "\n  [ MrBayes Partition Commands (Choose ONE) ]\n")
             # Generate MrBayes partition command(s)
             if gene_charset_names:
                 self.partition_text.insert(tk.END, "  partition by_gene = {0}: {1};\n".format(
                     len(gene_charset_names), " ".join(gene_charset_names))
                 )
             # Check if there are any codon names collected before generating the command
             if pos1_names or pos2_names or pos3_names:
                  partition_groups_mrbayes = []
                  if pos1_names: partition_groups_mrbayes.append(" ".join(pos1_names))
                  if pos2_names: partition_groups_mrbayes.append(" ".join(pos2_names))
                  if pos3_names: partition_groups_mrbayes.append(" ".join(pos3_names))
                  if partition_groups_mrbayes: # Ensure there's something to partition
                       self.partition_text.insert(tk.END, "  partition by_codon_pos = {0}: {1};\n".format(
                            len(partition_groups_mrbayes), ", ".join(partition_groups_mrbayes))
                       )
                  else:
                       self.partition_text.insert(tk.END, "  [ Cannot generate 'by_codon_pos' partition command - check DNA charsets ]\n")

             self.partition_text.insert(tk.END, "\n  [ Example: To use 'by_gene' partition: ]\n")
             if gene_charset_names:
                 self.partition_text.insert(tk.END, "  set partition = by_gene;\n")
             elif pos1_names or pos2_names or pos3_names: # Suggest codon partition if gene partition isn't possible but codon is
                  self.partition_text.insert(tk.END, "  [ Example: To use 'by_codon_pos' partition: ]\n")
                  self.partition_text.insert(tk.END, "  set partition = by_codon_pos;\n")
             else:
                 self.partition_text.insert(tk.END, "  [ No partitions generated ]\n")
             self.partition_text.insert(tk.END, "END; [mrbayes]\n")

        else:
             self.partition_text.insert(tk.END, "No partition data available from backend.\n")

        # Disable editing again
        self.concatenated_seqs_text.config(state="disabled")
        self.stats_text.config(state="disabled")
        self.partition_text.config(state="disabled")


    def clear_results_display(self):
        """Clears the text areas in the Results tab and stored results."""

        for text_widget in [self.concatenated_seqs_text, self.stats_text, self.partition_text]:
            text_widget.config(state="normal")
            text_widget.delete("1.0", tk.END)
            text_widget.config(state="disabled")
        self._concatenated_sequences = None
        self._partition_data = None
        self._statistics = None


    def on_reset(self):
        """Resets the input (genes list) and clears results."""
        # Only ask if there's something to reset (genes loaded or results present)
        if self.gene_names or self.raw_gene_contents or self._concatenated_sequences or self._partition_data or self._statistics:
            confirm = messagebox.askyesno("Confirm Reset", "Are you sure you want to remove all loaded genes and clear results?")
            if not confirm:
                return

        # Clear internal gene storage
        self.gene_names = []
        self.raw_gene_contents = []
        # Clear the listbox display
        self.gene_listbox.delete(0, tk.END)

        # Clear the results display and stored data
        self.clear_results_display()

    def export_fasta(self):
        """Exports the concatenated sequences to a FASTA file."""
        if not self._concatenated_sequences:
            messagebox.showwarning("Warning", "No concatenated sequences to export. Run 'Concatenate & Analyze' first.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=(("FASTA files", "*.fasta *.fas *.fna *.faa *.fa"), ("All files", "*.*")),
            title="Save Concatenated FASTA File"
        )

        if not filepath:
            return # User cancelled

        try:
            with open(filepath, "w") as f:
                # Sort by taxon name for consistent output
                sorted_taxons = sorted(self._concatenated_sequences.keys())
                for taxon in sorted_taxons:
                    sequence = self._concatenated_sequences[taxon]

                    # Add a newline at the end of the last line of sequence
                    wrapped_seq = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
                    f.write(f">{taxon}\n{wrapped_seq}\n") # Ensure a newline after the last sequence line

            messagebox.showinfo("Export Successful", f"Concatenated sequences exported to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save FASTA file:\n{e}")

    def export_partition(self):
        """Exports the partition data to a NEXUS formatted file (partition block only)."""
        if not self._partition_data:
            messagebox.showwarning("Warning", "No partition data to export. Run 'Concatenate & Analyze' first.")
            return
        # Check if we have sequence length info for the DATA block header
        num_taxa = len(self._concatenated_sequences) if self._concatenated_sequences else 0
        seq_len = len(list(self._concatenated_sequences.values())[0]) if self._concatenated_sequences else 0
        if num_taxa == 0 or seq_len == 0:
             messagebox.showwarning("Warning", "Concatenated sequence data is needed to determine matrix dimensions for the partition file header. Run 'Concatenate & Analyze' successfully first.")
             return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".nexus",
            filetypes=(("NEXUS files", "*.nexus *.nex"), ("Text files", "*.txt"), ("All files", "*.*")),
            title="Save Partition File (NEXUS partition block)"
        )

        if not filepath:
            return # User cancelled

        try:
            with open(filepath, "w") as f:
                # Generate NEXUS partition block similar to display_results, but without sequence data
                f.write("#nexus\n\n")


                f.write("BEGIN DATA;\n")
                # Add dimensions and format based on the generated data
                f.write(f"  DIMENSIONS NTAX={num_taxa} NCHAR={seq_len};\n")
                # Determine DATATYPE for the data block based on the partition data
                datatypes = set(item[2] for item in self._partition_data if item[2] != 'Unknown')
                if 'DNA' in datatypes and 'Protein' in datatypes:
                     # Mixed types, Nexus DATATYPE cannot represent this simply. Use DNA or UNKNOWN.
                     f.write("  FORMAT DATATYPE=DNA MISSING=- GAP=-;\n  [ WARNING: Mixed DNA/Protein data, DATATYPE=DNA may be unsuitable. ]\n")
                elif 'Protein' in datatypes:
                     f.write("  FORMAT DATATYPE=Protein MISSING=- GAP=-;\n")
                else: # Only DNA, only Unknown, or empty
                     f.write("  FORMAT DATATYPE=DNA MISSING=- GAP=-;\n  [ DATATYPE assumed as DNA ]\n")
                f.write("END; [DATA]\n\n")


                # ASSUMPTIONS or PAUP Block for CHARSETS and LINK
                # Using PAUP block as it supports both CHARSETS and LINK
                f.write("BEGIN PAUP;\n")
                f.write("  [ Charsets for gene partitions ]\n")
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                    clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                    f.write(f"  charset {clean_gene_name} = {coord_range};\n")

                # Add Codon Partitions for DNA genes
                f.write("\n  [ Charsets for codon positions (DNA genes only) ]\n")
                dna_genes_with_partitions_in_partition = [] # List of names of genes that got codon partitions
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                     if gene_type.upper() == 'DNA':
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          try:
                               start, end = map(int, coord_range.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    f.write(f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n")
                                    dna_genes_with_partitions_in_partition.append(clean_gene_name)
                               elif length > 0:
                                    f.write(f"  [ Warning: Gene {clean_gene_name} ({coord_range}) is too short ({length} bp) for full codon partitions. ]\n")
                          except ValueError:
                               f.write(f"  [ Warning: Could not parse range '{coord_range}' for codon partitions for gene {clean_gene_name}. ]\n")
                     elif gene_type.upper() == 'PROTEIN':
                          f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")


                if not dna_genes_with_partitions_in_partition:
                     f.write("  [ No DNA genes found or genes too short for codon partitions. ]\n")

                # Add LINK block
                f.write("\n  [ Link block (gene level) ]\n")
                if self._partition_data: # Only add LINK block if there's something to link
                     f.write("  link characters = ")
                     link_strings = []
                     for backend_gene_name, coord_range, gene_type in self._partition_data:
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          link_strings.append(f"{clean_gene_name} : {coord_range}")
                     f.write(", ".join(link_strings) + ";\n")
                else:
                     f.write("  [ No genes to link ]\n")



                if dna_genes_with_partitions_in_partition: # Only if DNA codon charsets were created
                    f.write("\n  [ Link block (by codon position, DNA genes only) ]\n")
                    f.write("  link characters = ")
                    codon_link_strings = []
                    # Recreate the range strings needed for the link block, filtering for genes that actually got codon partitions
                    pos1_ranges = [f"{re.sub(r'\W+', '_', name)}_pos1 : {start}-{end}\\3"
                                   for name, range_str, gene_type in self._partition_data if gene_type.upper() == 'DNA'
                                   for start, end in [map(int, range_str.split('-'))] if end-start+1 >= 3]
                    pos2_ranges = [f"{re.sub(r'\W+', '_', name)}_pos2 : {start+1}-{end}\\3"
                                   for name, range_str, gene_type in self._partition_data if gene_type.upper() == 'DNA'
                                   for start, end in [map(int, range_str.split('-'))] if end-start+1 >= 3]
                    pos3_ranges = [f"{re.sub(r'\W+', '_', name)}_pos3 : {start+2}-{end}\\3"
                                   for name, range_str, gene_type in self._partition_data if gene_type.upper() == 'DNA'
                                   for start, end in [map(int, range_str.split('-'))] if end-start+1 >= 3]

                    all_codon_ranges = pos1_ranges + pos2_ranges + pos3_ranges
                    if all_codon_ranges:
                        f.write(", ".join(all_codon_ranges) + ";\n")
                    else:
                         f.write("[ No DNA genes long enough for codon link block ]\n")


                f.write("END; [PAUP]\n\n")

                # Add MrBayes block example
                f.write("BEGIN mrbayes;\n")
                f.write("  [ Example MrBayes partition block ]\n")
                # Add charsets again for MrBayes block (both gene and codon levels)
                f.write("  [ Gene level charsets ]\n")
                gene_charset_names_mb = [] # Use a different list name to avoid confusion
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                    clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                    f.write(f"  charset {clean_gene_name} = {coord_range};\n")
                    gene_charset_names_mb.append(clean_gene_name)


                f.write("\n  [ Codon position charsets (DNA genes only) ]\n")
                mr_bayes_dna_charsets_exist_mb = False
                pos1_names_mb = []
                pos2_names_mb = []
                pos3_names_mb = []

                for backend_gene_name, coord_range, gene_type in self._partition_data:
                     if gene_type.upper() == 'DNA':
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          try:
                               start, end = map(int, coord_range.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    f.write(f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n")
                                    mr_bayes_dna_charsets_exist_mb = True
                                    pos1_names_mb.append(f"{clean_gene_name}_pos1")
                                    pos2_names_mb.append(f"{clean_gene_name}_pos2")
                                    pos3_names_mb.append(f"{clean_gene_name}_pos3")
                               # Warning about short genes is handled in PAUP block display
                          except ValueError:
                               pass # Error message handled in PAUP block display
                     elif gene_type.upper() == 'PROTEIN':
                          f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")


                if not mr_bayes_dna_charsets_exist_mb:
                     f.write("  [ No DNA genes found or genes too short for codon partitions. ]\n")


                f.write("\n  [ MrBayes Partition Commands (Choose ONE) ]\n")
                # Generate MrBayes partition command(s)
                if gene_charset_names_mb:
                    f.write("  partition by_gene = {0}: {1};\n".format(
                        len(gene_charset_names_mb), " ".join(gene_charset_names_mb))
                    )
                # Check if any codon names were collected before generating the command
                if pos1_names_mb or pos2_names_mb or pos3_names_mb:
                     partition_groups_mrbayes_mb = []
                     if pos1_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos1_names_mb))
                     if pos2_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos2_names_mb))
                     if pos3_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos3_names_mb))
                     if partition_groups_mrbayes_mb: # Ensure there's something to partition
                          f.write("  partition by_codon_pos = {0}: {1};\n".format(
                               len(partition_groups_mrbayes_mb), ", ".join(partition_groups_mrbayes_mb))
                          )
                     else:
                          f.write("  [ Cannot generate 'by_codon_pos' partition command - check DNA charsets ]\n")

                f.write("\n  [ Example: To use 'by_gene' partition: ]\n")
                if gene_charset_names_mb:
                    f.write("  set partition = by_gene;\n")
                elif pos1_names_mb or pos2_names_mb or pos3_names_mb: # Suggest codon partition if gene partition isn't possible but codon is
                     f.write("  [ Example: To use 'by_codon_pos' partition: ]\n")
                     f.write("  set partition = by_codon_pos;\n")
                else:
                    f.write("  [ No partitions generated ]\n")

                f.write("END; [mrbayes]\n")


            messagebox.showinfo("Export Successful", f"Partition data exported to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save partition file:\n{e}")

    def export_full_nexus(self):
        """Exports the concatenated sequences and partition data to a full NEXUS file."""
        if not self._concatenated_sequences:
            messagebox.showwarning("Warning", "No concatenated sequences to export. Run 'Concatenate & Analyze' first.")
            return
        if not self._partition_data:
             # Partition data is needed for CHARSETS/LINK blocks
             messagebox.showwarning("Warning", "No partition data available. Run 'Concatenate & Analyze' successfully first.")
             return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".nexus",
            filetypes=(("NEXUS files", "*.nexus *.nex"), ("All files", "*.*")),
            title="Save Full NEXUS File"
        )

        if not filepath:
            return # User cancelled

        try:
            with open(filepath, "w") as f:
                f.write("#nexus\n\n")

                num_taxa = len(self._concatenated_sequences)
                seq_len = len(list(self._concatenated_sequences.values())[0]) # Assumes all concatenated sequences are same length

                f.write("BEGIN TAXA;\n")
                f.write(f"  DIMENSIONS NTAX={num_taxa};\n")
                f.write("  TAXLABELS\n")
                sorted_taxons = sorted(self._concatenated_sequences.keys())
                for taxon in sorted_taxons:
                    # Add quotes around taxon names if they contain spaces or special characters
                    # Check for standard problematic chars plus whitespace
                    if re.search(r'[\s\'"`=;:,]+', str(taxon)):
                         f.write(f"    '{taxon}'\n")
                    else:
                         f.write(f"    {taxon}\n")
                f.write("  ;\n") # End TAXLABELS
                f.write("END; [TAXA]\n\n")


                # DATA Block
                f.write("BEGIN DATA;\n")
                f.write(f"  DIMENSIONS NTAX={num_taxa} NCHAR={seq_len};\n")

                # Determine DATATYPE for the data block based on the partition data
                datatypes = set(item[2] for item in self._partition_data if item[2] != 'Unknown')
                if 'DNA' in datatypes and 'Protein' in datatypes:
                     # Mixed types, Nexus DATATYPE cannot represent this simply. Use DNA or UNKNOWN.
                     f.write("  FORMAT DATATYPE=DNA MISSING=- GAP=- INTERLEAVE=NO;\n  [ WARNING: Mixed DNA/Protein data, DATATYPE=DNA may be unsuitable. ]\n")
                elif 'Protein' in datatypes:
                     f.write("  FORMAT DATATYPE=Protein MISSING=- GAP=- INTERLEAVE=NO;\n")
                else: # Only DNA, only Unknown, or empty
                     f.write("  FORMAT DATATYPE=DNA MISSING=- GAP=- INTERLEAVE=NO;\n  [ DATATYPE assumed as DNA ]\n")


                f.write("  MATRIX;\n")

                # Write sequences (sequential format)
                for taxon in sorted_taxons:
                    sequence = self._concatenated_sequences[taxon]
                    # Add quotes around taxon names if they contain spaces or special characters for NEXUS
                    taxon_label = f"'{taxon}'" if re.search(r'[\s\'"`=;:,]+', str(taxon)) else taxon
                    # Write label and sequence
                    f.write(f"    {taxon_label}\n")
                    # Write sequence, possibly wrapped
                    wrapped_seq = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)]) # Wrap at 60 for readability
                    # Indent wrapped lines
                    indented_wrapped_seq = '\n'.join([f"      {line}" for line in wrapped_seq.splitlines()])
                    f.write(f"{indented_wrapped_seq}\n")

                f.write("  ;\n") # End MATRIX
                f.write("END; [DATA]\n\n")


                # ASSUMPTIONS or PAUP Block for CHARSETS and LINK
                # Using PAUP block as it supports both CHARSETS and LINK
                f.write("BEGIN PAUP;\n")
                f.write("  [ Charsets for gene partitions ]\n")
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                    clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                    f.write(f"  charset {clean_gene_name} = {coord_range};\n")

                # Add Codon Partitions for DNA genes
                f.write("\n  [ Charsets for codon positions (DNA genes only) ]\n")
                dna_genes_with_partitions_in_nexus = [] # Track which DNA genes got codon partitions
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                     if gene_type.upper() == 'DNA':
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          try:
                               start, end = map(int, coord_range.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    f.write(f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n")
                                    dna_genes_with_partitions_in_nexus.append(clean_gene_name)
                               elif length > 0:
                                    f.write(f"  [ Warning: Gene {clean_gene_name} ({coord_range}) is too short ({length} bp) for full codon partitions. ]\n")
                          except ValueError:
                               f.write(f"  [ Warning: Could not parse range '{coord_range}' for codon partitions for gene {clean_gene_name}. ]\n")
                     elif gene_type.upper() == 'PROTEIN':
                           f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")


                if not dna_genes_with_partitions_in_nexus:
                     f.write("  [ No DNA genes found or genes too short for codon partitions. ]\n")


                # Add LINK block
                f.write("\n  [ Link block (gene level) ]\n")
                if self._partition_data: # Only add LINK block if there's something to link
                    f.write("  link characters = ")
                    link_strings = []
                    for backend_gene_name, coord_range, gene_type in self._partition_data:
                         clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                         link_strings.append(f"{clean_gene_name} : {coord_range}")
                    f.write(", ".join(link_strings) + ";\n")
                else:
                    f.write("  [ No genes to link ]\n")


                if dna_genes_with_partitions_in_nexus: # Only if DNA codon charsets were created
                    f.write("\n  [ Link block (by codon position, DNA genes only) ]\n")
                    f.write("  link characters = ")

                     # Recreate the range strings needed for the link block, filtering for genes that actually got codon partitions
                    pos1_ranges = [f"{re.sub(r'\W+', '_', name)}_pos1 : {start}-{end}\\3"
                                   for name, range_str, gene_type in self._partition_data if gene_type.upper() == 'DNA'
                                   for start, end in [map(int, range_str.split('-'))] if end-start+1 >= 3]
                    pos2_ranges = [f"{re.sub(r'\W+', '_', name)}_pos2 : {start+1}-{end}\\3"
                                   for name, range_str, gene_type in self._partition_data if gene_type.upper() == 'DNA'
                                   for start, end in [map(int, range_str.split('-'))] if end-start+1 >= 3] # Corrected typo: range_str instead of range_range
                    pos3_ranges = [f"{re.sub(r'\W+', '_', name)}_pos3 : {start+2}-{end}\\3"
                                   for name, range_str, gene_type in self._partition_data if gene_type.upper() == 'DNA'
                                   for start, end in [map(int, range_str.split('-'))] if end-start+1 >= 3]


                    all_codon_ranges = pos1_ranges + pos2_ranges + pos3_ranges
                    if all_codon_ranges:
                        f.write(", ".join(all_codon_ranges) + ";\n")
                    else:
                         f.write("[ No DNA genes long enough for codon link block ]\n")

                f.write("END; [PAUP]\n\n")

                # Add MrBayes block example
                f.write("BEGIN mrbayes;\n")
                f.write("  [ Example MrBayes partition block ]\n")
                # Add charsets again for MrBayes block (both gene and codon levels)
                f.write("  [ Gene level charsets ]\n")
                gene_charset_names_mb = [] # Use a different list name to avoid confusion
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                    clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                    f.write(f"  charset {clean_gene_name} = {coord_range};\n")
                    gene_charset_names_mb.append(clean_gene_name)


                f.write("\n  [ Codon position charsets (DNA genes only) ]\n")
                mr_bayes_dna_charsets_exist_mb = False
                pos1_names_mb = []
                pos2_names_mb = []
                pos3_names_mb = []

                for backend_gene_name, coord_range, gene_type in self._partition_data:
                     if gene_type.upper() == 'DNA':
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          try:
                               start, end = map(int, coord_range.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    f.write(f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n")
                                    mr_bayes_dna_charsets_exist_mb = True
                                    pos1_names_mb.append(f"{clean_gene_name}_pos1")
                                    pos2_names_mb.append(f"{clean_gene_name}_pos2")
                                    pos3_names_mb.append(f"{clean_gene_name}_pos3")
                               # Warning about short genes is handled in PAUP block display
                          except ValueError:
                               pass # Error message handled in PAUP block display
                     elif gene_type.upper() == 'PROTEIN':
                          f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")


                if not mr_bayes_dna_charsets_exist_mb:
                     f.write("  [ No DNA genes found or genes too short for codon partitions. ]\n")


                f.write("\n  [ MrBayes Partition Commands (Choose ONE) ]\n")
                # Generate MrBayes partition command(s)
                if gene_charset_names_mb:
                    f.write("  partition by_gene = {0}: {1};\n".format(
                        len(gene_charset_names_mb), " ".join(gene_charset_names_mb))
                    )
                # Check if any codon names were collected before generating the command
                if pos1_names_mb or pos2_names_mb or pos3_names_mb:
                     partition_groups_mrbayes_mb = []
                     if pos1_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos1_names_mb))
                     if pos2_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos2_names_mb))
                     if pos3_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos3_names_mb))
                     if partition_groups_mrbayes_mb: # Ensure there's something to partition
                          f.write("  partition by_codon_pos = {0}: {1};\n".format(
                               len(partition_groups_mrbayes_mb), ", ".join(partition_groups_mrbayes_mb))
                          )
                     else:
                          f.write("  [ Cannot generate 'by_codon_pos' partition command - check DNA charsets ]\n")

                f.write("\n  [ Example: To use 'by_gene' partition: ]\n")
                if gene_charset_names_mb:
                    f.write("  set partition = by_gene;\n")
                elif pos1_names_mb or pos2_names_mb or pos3_names_mb: # Suggest codon partition if gene partition isn't possible but codon is
                     f.write("  [ Example: To use 'by_codon_pos' partition: ]\n")
                     f.write("  set partition = by_codon_pos;\n")
                else:
                    f.write("  [ No partitions generated ]\n")

                f.write("END; [mrbayes]\n")


            messagebox.showinfo("Export Successful", f"Full NEXUS file exported to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save Full NEXUS file:\n{e}")

# --- Main execution ---
if __name__ == "__main__":

    root = tk.Tk()
    app = SequenceConcatenatorApp(root)
    root.mainloop()