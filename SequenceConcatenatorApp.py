import sys
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import re
import os
import traceback
from turtledemo.penrose import f

# Assuming SequenceConcatenator.py is in the same directory or accessible in the python path
try:
    # Pass the assigned gene names along with raw data to the backend
    from SequenceConcatenator import SequenceConcatenator
except ImportError:
    messagebox.showerror("Import Error", "Could not find SequenceConcatenator.py. Please ensure it is in the same directory.")
    exit()


class SequenceConcatenatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Sequence Concatenator UI")
        self.root.geometry("1200x900")
        self.root.configure(bg="#f4f4f4")

        self.gene_names = [] # Stores the names assigned by the frontend
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
        self.font_normal = ("Arial", 10)
        self.font_bold = ("Arial", 10, "bold")
        self.font_heading = ("Arial", 11, "bold")
        self.font_button = ("Segoe UI", 10, "bold")
        self.font_mono = ("Courier New", 9)


        # --- Menu Bar (Removed) ---
        # The menubar creation code has been removed as requested previously.


        # --- Toolbar (Removed) ---
        # The toolbar frame and its buttons have been removed as requested previously.


        # Notebook (Tabs)
        self.notebook = ttk.Notebook(root)
        # Added padding=0 to remove default padding at the top, effectively moving content up
        self.notebook.pack(pady=0, padx=10, fill="both", expand=True)

        # Create tabs
        self.tabs = {
            "Input": ttk.Frame(self.notebook, padding="10", style='TFrame'),
            "Results": ttk.Frame(self.notebook, padding="10", style='TFrame'),
        }

        # Configure ttk Style
        style = ttk.Style()
        style.theme_use('default')
        style.configure('TFrame', background=self.bg_color)
        style.configure('TNotebook.Tab', padding=[10, 5], font=self.font_normal)
        style.map('TNotebook.Tab',
                  background=[('selected', self.primary_color)],
                  foreground=[('selected', 'white')],
                  expand=[('selected', [1, 1, 1, 0])])

        # Configure Treeview style
        style.configure("Treeview.Heading", font=self.font_bold, background="#eeeeee", foreground="black")
        style.configure("Treeview", font=self.font_normal, rowheight=20)
        style.map("Treeview", background=[("selected", self.primary_color)], foreground=[("selected", "white")],
                 font=[('selected', self.font_normal)]) # Keep normal font on select unless editing

        # Custom tag for Outgroup row - Apply bold font and different background
        # Use a different background color when *not* selected to make it stand out
        style.configure("Treeview.outgroup", font=self.font_bold)
        style.map("Treeview.outgroup",
                  background=[("selected", self.primary_color), ('!selected', '#d9d9d9')], # Light gray background when not selected
                  foreground=[("selected", "white"), ('!selected', "black")], # Black text when not selected
                  font=[('selected', self.font_bold), ('!selected', self.font_bold)]) # Always bold font for outgroup row


        for tab_name, tab_frame in self.tabs.items():
            self.notebook.add(tab_frame, text=tab_name)

        self.setup_input_tab()

        # Container frame for results tab content
        self.results_content_frame = ttk.Frame(self.tabs["Results"], style='TFrame')
        self.results_content_frame.pack(fill="both", expand=True)
        self.setup_results_tab(self.results_content_frame)


        # Store results internally for export and display
        self._concatenated_sequences = None # {taxon_name: sequence_string, ...} -> Keys are potentially edited
        # Store the SequenceConcatenator instance
        self._concatenator = None
        # Store gene info processed by the backend (needed for recalculating divergence)
        self._backend_gene_info = None
        self._partition_data = None # [(gene_name, 'start-end', gene_type), ...]
        self._statistics = None # Dictionary holding various stats, including divergence (initial)
        self._divergence_data = None # {taxon_name: {'Total score': float, 'No of charsets': int, 'ActualGeneName': 'count (perc%)', ...}, ...}

        # Variable to hold the current reference taxon for divergence calculation
        self._reference_taxon = None

        # Variable to hold the Treeview Entry widget for editing (Keeping taxon renaming)
        self.entry_editor = None
        self.entry_item = None
        self.entry_column = None
        # Map to track Treeview item IDs to current taxon names for editing/deleting (Keep editing)
        self._treeview_id_to_taxon_name = {}


    def setup_input_tab(self):
        input_frame = self.tabs["Input"]
        input_frame.columnconfigure(0, weight=1)
        input_frame.rowconfigure(0, weight=1)

        list_container_frame = tk.Frame(input_frame, bg=self.bg_color)
        list_container_frame.grid(row=0, column=0, pady=10, padx=10, sticky="nsew")

        tk.Label(list_container_frame, text="Loaded Genes:", font=self.font_heading, bg=self.bg_color, fg=self.heading_color).pack(anchor="w", pady=(0, 5))

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
            selectmode=tk.EXTENDED
        )
        self.gene_listbox.pack(side="left", fill="both", expand=True)

        list_scrollbar = tk.Scrollbar(listbox_frame, orient="vertical", command=self.gene_listbox.yview)
        list_scrollbar.pack(side="right", fill="y")
        self.gene_listbox.config(yscrollcommand=list_scrollbar.set)

        btn_frame = tk.Frame(input_frame, bg=self.bg_color)
        btn_frame.grid(row=1, column=0, sticky="ew", pady=10, padx=10)

        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)
        btn_frame.columnconfigure(2, weight=0)
        btn_frame.columnconfigure(3, weight=0)

        # Kept Add/Remove/Reset/Submit buttons on Input tab
        add_gene_btn = tk.Button(
            btn_frame,
            text="Add Gene (Manual/File)",
            command=self.open_add_gene_dialog,
            bg=self.button_bg_primary, fg=self.button_fg,
            activebackground=self.active_bg_primary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        add_gene_btn.grid(row=0, column=0, sticky="w", padx=(0, 10))

        remove_gene_btn = tk.Button(
            btn_frame,
            text="Remove Selected Gene(s)",
            command=self.remove_selected_gene,
            bg=self.button_bg_secondary, fg=self.button_fg,
            activebackground=self.active_bg_secondary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        remove_gene_btn.grid(row=0, column=1, sticky="w", padx=(0, 10))

        action_btn_frame = tk.Frame(btn_frame, bg=self.bg_color)
        action_btn_frame.grid(row=0, column=3, sticky="e")

        reset_btn = tk.Button(
            action_btn_frame,
            text="Reset All",
            command=self.on_reset,
            bg=self.button_bg_secondary, fg=self.button_fg,
            activebackground=self.active_bg_secondary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        reset_btn.pack(side="left", padx=(0, 10))

        submit_btn = tk.Button(
            action_btn_frame,
            text="Concatenate & Analyze",
            command=self.on_submit,
            bg=self.button_bg_primary, fg=self.button_fg,
            activebackground=self.active_bg_primary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        submit_btn.pack(side="left")


    def setup_results_tab(self, parent_frame):
        parent_frame.columnconfigure(0, weight=1)
        parent_frame.rowconfigure(0, weight=0) # Top panes (stats/partition) fixed height
        parent_frame.rowconfigure(1, weight=1) # Divergence table takes main space
        parent_frame.rowconfigure(2, weight=0) # Buttons fixed height


        # Frame for Statistics and Partition
        top_panes_frame = ttk.Frame(parent_frame, style='TFrame')
        top_panes_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        top_panes_frame.columnconfigure(0, weight=1)
        top_panes_frame.columnconfigure(1, weight=1)
        top_panes_frame.rowconfigure(0, weight=1)


        # Frame for Statistics Display ---
        stats_frame = tk.LabelFrame(top_panes_frame, text="Statistics", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=5)
        stats_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        stats_frame.columnconfigure(0, weight=1)
        stats_frame.rowconfigure(0, weight=1)

        self.stats_text = tk.Text(
            stats_frame, font=self.font_normal, bg="white", fg=self.text_color,
            height=8, state="disabled", wrap="word"
        )
        self.stats_text.grid(row=0, column=0, sticky="nsew")

        # Frame for Partition Display ---
        partition_frame = tk.LabelFrame(top_panes_frame, text="Partition (NEXUS format example)", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=5)
        partition_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)
        partition_frame.columnconfigure(0, weight=1)
        partition_frame.rowconfigure(0, weight=1)

        partition_text_frame = tk.Frame(partition_frame)
        partition_text_frame.grid(row=0, column=0, sticky="nsew")
        partition_text_frame.columnconfigure(0, weight=1)
        partition_text_frame.rowconfigure(0, weight=1)

        self.partition_text = tk.Text(
            partition_text_frame, font=self.font_mono, bg="white", fg=self.text_color,
            height=8, state="disabled", wrap="none"
        )
        self.partition_text.grid(row=0, column=0, sticky="nsew")

        partition_scrollbar_y = tk.Scrollbar(partition_text_frame, orient="vertical", command=self.partition_text.yview)
        partition_scrollbar_y.grid(row=0, column=1, sticky="ns")
        self.partition_text.config(yscrollcommand=partition_scrollbar_y.set)

        partition_scrollbar_x = tk.Scrollbar(partition_frame, orient="horizontal", command=self.partition_text.xview)
        partition_scrollbar_x.grid(row=1, column=0, sticky="ew")
        self.partition_text.config(xscrollcommand=partition_scrollbar_x.set)


        # --- Divergence Table ---
        divergence_frame = tk.LabelFrame(parent_frame, text="Divergence from Reference Taxon", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=5)
        divergence_frame.grid(row=1, column=0, sticky="nsew", padx=5, pady=5)
        divergence_frame.columnconfigure(0, weight=1)
        divergence_frame.rowconfigure(0, weight=1)


        # Frame for Treeview and Scrollbars
        treeview_frame = tk.Frame(divergence_frame)
        treeview_frame.grid(row=0, column=0, sticky="nsew")
        treeview_frame.columnconfigure(0, weight=1)
        treeview_frame.rowconfigure(0, weight=1)


        self.divergence_treeview = ttk.Treeview(treeview_frame, show="headings", selectmode="browse")
        self.divergence_treeview.grid(row=0, column=0, sticky="nsew")

        tree_scrollbar_y = tk.Scrollbar(treeview_frame, orient="vertical", command=self.divergence_treeview.yview)
        tree_scrollbar_y.grid(row=0, column=1, sticky="ns")
        self.divergence_treeview.config(yscrollcommand=tree_scrollbar_y.set)

        tree_scrollbar_x = tk.Scrollbar(divergence_frame, orient="horizontal", command=self.divergence_treeview.xview)
        tree_scrollbar_x.grid(row=1, column=0, sticky="ew")
        self.divergence_treeview.config(xscrollcommand=tree_scrollbar_x.set)


        # Bind click for editing (Keeping taxon renaming)
        self.divergence_treeview.bind('<ButtonRelease-1>', self._on_treeview_click)


        # --- Buttons at bottom of Results tab ---
        button_frame = tk.Frame(parent_frame, bg=self.bg_color)
        button_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10), padx=10)

        # Re-added Export buttons
        export_fasta_btn = tk.Button(
            button_frame, text="Export FASTA", command=self.export_fasta,
            bg=self.button_bg_primary, fg=self.button_fg,
            activebackground=self.active_bg_primary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_fasta_btn.pack(side="left", padx=(0, 10))

        export_partition_btn = tk.Button(
            button_frame, text="Export Partition (NEXUS block)", command=self.export_partition,
            bg=self.button_bg_primary, fg=self.button_fg,
            activebackground=self.active_bg_primary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_partition_btn.pack(side="left", padx=(0, 10))

        export_full_nexus_btn = tk.Button(
             button_frame, text="Export Full NEXUS", command=self.export_full_nexus,
             bg=self.button_bg_primary, fg=self.button_fg,
             activebackground=self.active_bg_primary, activeforeground=self.button_fg,
             font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_full_nexus_btn.pack(side="left") # No padx on the last button on the left

        # Add the "Make Outgroup" button - placed on the right
        make_outgroup_btn = tk.Button(
            button_frame, text="Make Outgroup", command=self.make_selected_taxon_outgroup,
            bg=self.button_bg_secondary, fg=self.button_fg, # Using secondary color for this action button
            activebackground=self.active_bg_secondary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        make_outgroup_btn.pack(side="right") # Place on the right side


    # Kept as this is the primary way to add genes
    def open_add_gene_dialog(self):
        """Opens a new window to add a gene manually or load from file."""
        dialog = tk.Toplevel(self.root)
        dialog.title("Add New Gene or Load File")
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.configure(bg=self.bg_color)
        main_x = self.root.winfo_x()
        main_y = self.root.winfo_y()
        main_width = self.root.winfo_width()
        main_height = self.root.winfo_height()
        dialog_width = 600
        dialog_height = 700
        dialog_x = main_x + (main_width // 2) - (dialog_width // 2)
        dialog_y = main_y + (main_height // 2) - (dialog_height // 2)
        dialog.geometry(f"{dialog_width}x{dialog_height}+{dialog_x}+{dialog_y}")
        dialog.resizable(True, True)

        dialog_frame = tk.Frame(dialog, bg=self.bg_color)
        dialog_frame.pack(fill="both", expand=True, padx=15, pady=10)
        dialog_frame.columnconfigure(0, weight=1)
        dialog_frame.rowconfigure(0, weight=3)
        dialog_frame.rowconfigure(1, weight=0)
        dialog_frame.rowconfigure(2, weight=0)

        manual_frame = tk.LabelFrame(dialog_frame, text="Add Manually (Paste Multi-FASTA Data)", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=10)
        manual_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 10))
        manual_frame.columnconfigure(0, weight=1)
        manual_frame.rowconfigure(0, weight=0)
        manual_frame.rowconfigure(1, weight=0)
        manual_frame.rowconfigure(2, weight=0)
        manual_frame.rowconfigure(3, weight=1)
        manual_frame.rowconfigure(4, weight=0)

        tk.Label(manual_frame, text="Gene Name :", font=self.font_normal, bg=self.bg_color, fg=self.text_color).grid(row=0, column=0, sticky="w")
        name_entry = tk.Entry(manual_frame, font=self.font_normal)
        name_entry.grid(row=1, column=0, sticky="ew", pady=(0, 5))

        tk.Label(manual_frame, text="Sequence Data (Paste Multi-FASTA here):", font=self.font_normal, bg=self.bg_color, fg=self.text_color).grid(row=2, column=0, sticky="nw")

        manual_text_frame = tk.Frame(manual_frame)
        manual_text_frame.grid(row=3, column=0, sticky="nsew", pady=(0, 5))
        manual_text_frame.columnconfigure(0, weight=1)
        manual_text_frame.rowconfigure(0, weight=1)

        seq_text = tk.Text(manual_text_frame, font=self.font_mono, wrap="word")
        seq_text.grid(row=0, column=0, sticky="nsew")
        manual_scrollbar_y = tk.Scrollbar(manual_text_frame, orient="vertical", command=seq_text.yview)
        manual_scrollbar_y.grid(row=0, column=1, sticky="ns")
        seq_text.config(yscrollcommand=manual_scrollbar_y.set)

        add_manual_btn = tk.Button(
            manual_frame, text="Add Manual Gene",
            command=lambda: self.add_gene_from_dialog(dialog, name_entry.get(), seq_text.get("1.0", "end").strip()),
            bg=self.button_bg_primary, fg=self.button_fg,
            activebackground=self.active_bg_primary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=10, pady=5, cursor="hand2"
        )
        add_manual_btn.grid(row=4, column=0, pady=(5,0))


        load_frame = tk.LabelFrame(dialog_frame, text="Load from File(s) (FASTA, NEXUS, GenBank)", font=self.font_heading, bg=self.bg_color, fg=self.heading_color, padx=10, pady=10)
        load_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        load_frame.columnconfigure(0, weight=1)
        load_frame.rowconfigure(0, weight=0)
        load_frame.rowconfigure(1, weight=0)
        load_frame.rowconfigure(2, weight=0)

        tk.Label(load_frame, text="Gene Name for File(s) (Optional - uses filename if empty):", font=self.font_normal, bg=self.bg_color, fg=self.text_color).grid(row=0, column=0, sticky="w")
        load_name_entry = tk.Entry(load_frame, font=self.font_normal)
        load_name_entry.grid(row=1, column=0, sticky="ew", pady=(0, 5))

        load_file_btn = tk.Button(
            load_frame, text="Browse and Load File(s)",
            command=lambda: self.load_file_dialog_and_store_raw(dialog, load_name_entry.get().strip()),
            bg=self.button_bg_primary, fg=self.button_fg,
            activebackground=self.active_bg_primary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=10, pady=5, cursor="hand2"
        )
        load_file_btn.grid(row=2, column=0, pady=(5,0))


        close_btn = tk.Button(
            dialog_frame, text="Cancel", command=dialog.destroy,
            bg=self.button_bg_secondary, fg=self.button_fg,
            activebackground=self.active_bg_secondary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=10, pady=5, cursor="hand2"
        )
        close_btn.grid(row=2, column=0, pady=10)

        dialog.focus_set()

    # Kept as it's used by the add gene dialog
    def add_gene_from_dialog(self, dialog, gene_name, raw_content):
        gene_name = gene_name.strip()
        raw_content = raw_content.strip()

        if not gene_name:
            messagebox.showwarning("Warning", "Please provide a gene name.", parent=dialog)
            return
        if not raw_content:
            messagebox.showwarning("Warning", "Please paste sequence data.", parent=dialog)
            return

        raw_lines = raw_content.splitlines()
        self._add_raw_gene_data(gene_name, raw_lines)
        dialog.destroy()

    # Kept as it's used by the add gene dialog
    def load_file_dialog_and_store_raw(self, dialog, manual_gene_name=""):
        filepaths = filedialog.askopenfilenames(
            initialdir=".",
            title="Select Gene File(s)",
            filetypes=(
                ("Sequence files", "*.fasta *.fas *.fna *.faa *.fa *.nexus *.nex *.gbff"),
                ("FASTA files", "*.fasta *.fas *.fna *.faa *.fa"),
                ("NEXUS files", "*.nexus *.nex"),
                ("GenBank files", "*.gbff"),
                ("All files", "*.*")
            ),
            parent=dialog
        )
        if not filepaths:
            return

        successfully_loaded_names = []
        failed_files = []

        for filepath in filepaths:
            try:
                with open(filepath, "r", encoding='utf-8', errors='ignore') as f:
                    raw_content_lines = f.readlines()

                if not raw_content_lines or all(line.strip() == '' for line in raw_content_lines):
                    messagebox.showwarning("Warning", f"File '{os.path.basename(filepath)}' is empty or contains only whitespace. Skipping.", parent=dialog)
                    failed_files.append(f"{os.path.basename(filepath)} (Empty)")
                    continue

                gene_name = manual_gene_name if manual_gene_name else os.path.splitext(os.path.basename(filepath))[0]

                added_name = self._add_raw_gene_data(gene_name, raw_content_lines)
                successfully_loaded_names.append(added_name)

            except Exception as e:
                messagebox.showerror("File Reading Error", f"Error reading file '{os.path.basename(filepath)}':\n{e}", parent=dialog)
                failed_files.append(f"{os.path.basename(filepath)} (Error: {e})")


        status_msg = ""
        if successfully_loaded_names:
             status_msg += "Successfully loaded the following gene file(s) (parsing happens on submit):\n" + "\n".join(successfully_loaded_names)
        if failed_files:
             if status_msg: status_msg += "\n\n"
             status_msg += "Failed to load the following file(s):\n" + "\n".join(failed_files)

        if status_msg:
             if successfully_loaded_names:
                 messagebox.showinfo("Load Status", status_msg, parent=dialog)
                 dialog.destroy() # Close dialog only on successful load
             else:
                  messagebox.showerror("Load Failed", status_msg, parent=dialog)
        else:
             # No files selected, no files loaded, no files failed
             pass


    # Kept as it handles adding to the internal lists and listbox
    def _add_raw_gene_data(self, gene_name, raw_lines):
        if not gene_name or not raw_lines:
            return

        original_name = gene_name
        counter = 1
        existing_names_in_listbox = list(self.gene_listbox.get(0, tk.END))
        temp_name = gene_name
        while temp_name in existing_names_in_listbox:
            temp_name = f"{original_name}_{counter}"
            counter += 1
        gene_name = temp_name

        self.gene_names.append(gene_name)
        self.raw_gene_contents.append(raw_lines)
        self.gene_listbox.insert(tk.END, gene_name)

        return gene_name

    # Kept as it's triggered by the "Remove Selected Gene(s)" button on the Input tab
    def remove_selected_gene(self):
        selected_indices = self.gene_listbox.curselection()
        if not selected_indices:
            messagebox.showwarning("Warning", "Please select one or more genes to remove.")
            return

        selected_gene_names = [self.gene_listbox.get(i) for i in selected_indices]
        confirm = messagebox.askyesno("Confirm Removal", f"Are you sure you want to remove the selected gene(s)?\n\n" + "\n".join(selected_gene_names))
        if not confirm:
            return

        # Remove from lists and listbox in reverse order to maintain correctness
        for i in sorted(selected_indices, reverse=True):
            self.gene_listbox.delete(i)
            del self.gene_names[i]
            del self.raw_gene_contents[i]

        self.clear_results_display()

    # Kept as it's triggered by the "Concatenate & Analyze" button on the Input tab
    def on_submit(self):
        if not self.raw_gene_contents:
            messagebox.showwarning("Warning", "Please add at least one gene before submitting.")
            return
        # Ensure gene_names and raw_gene_contents lists are in sync
        if len(self.gene_names) != len(self.raw_gene_contents):
             messagebox.showerror("Internal Error", "Mismatch between gene names and raw contents list lengths.")
             return


        self.clear_results_display()
        self._treeview_id_to_taxon_name = {} # Clear map for new results
        self._reference_taxon = None # Reset reference taxon on new analysis
        self._concatenator = None # Clear previous backend instance
        self._backend_gene_info = None # Clear previous backend gene info


        try:
            # Create and store the SequenceConcatenator instance
            self._concatenator = SequenceConcatenator(self.raw_gene_contents, self.gene_names)

            # Get results from the instance
            self._concatenated_sequences = self._concatenator.get_concatenated_sequences()
            self._backend_gene_info = self._concatenator.get_processed_gene_info() # Get processed gene info
            self._partition_data = self._concatenator.get_partition()
            # Statistics are calculated by the backend using its default reference (first taxon)
            self._statistics = self._concatenator.get_statistics()


            # Initial divergence data comes from the initial statistics calculation
            self._divergence_data = self._statistics.get("Divergence Data", {})

            if not self._concatenated_sequences:
                 messagebox.showwarning("Concatenation Warning", "No concatenated sequences were produced by the backend. This might happen if no valid data was parsed or there are no common taxa across genes.")
                 self.clear_results_display()
                 return

            # Set the initial reference taxon to the one used by the backend (the first one in the sorted concatenated keys)
            taxons_in_concat = sorted(list(self._concatenated_sequences.keys()))
            if taxons_in_concat:
                 # The backend calculated initial divergence based on the first taxon in its sorted list
                 self._reference_taxon = taxons_in_concat[0]
            else:
                 self._reference_taxon = None # Should not happen if concat_sequences is not empty


            # Check if divergence data seems incomplete when taxons exist
            if len(taxons_in_concat) > 1 and not self._divergence_data:
                 messagebox.showwarning("Divergence Warning", "No divergence data was produced by the backend. Ensure your data has sufficient overlaps between taxons and genes.")
            elif len(taxons_in_concat) > 1 and self._divergence_data and len(self._divergence_data) <= 1:
                 # This case happens if only the reference taxon's row exists in divergence data
                 messagebox.showwarning("Divergence Warning", f"Divergence data calculated only for reference taxon '{self._reference_taxon}'. Ensure other taxons have sufficient data.")


            # Display the results, including the initial divergence table
            self.display_results_panels(self._statistics, self._partition_data)
            # Pass the initial divergence data calculated by the backend
            self._update_divergence_table_display(self._divergence_data, self._reference_taxon)


            self.notebook.select(self.tabs["Results"])

        except Exception as e:
            messagebox.showerror("Processing Error", f"An error occurred during sequence processing.\nError details: {e}")
            traceback.print_exc()
            self.clear_results_display()

    # New method to handle displaying results panels (Stats and Partition)
    def display_results_panels(self, statistics, partition):
        """Populates the Statistics and Partition text areas."""
         # --- Display Statistics ---
        self.stats_text.config(state="normal")
        self.stats_text.delete("1.0", tk.END)
        self.stats_text.insert(tk.END, "Statistics:\n\n")
        if statistics:
            preferred_order = [
                "Number of Taxa",
                "Number of Genes",
                "Total Length",
                "Percentage Overall Missing Data (%)",
                "Taxon-Gene Matrix Sparsity (%)",
                "Missing Data per Taxon (count)",
                "Missing Data per Gene (count)",
                "Gene Lengths (in concatenated alignment)",
            ]
            displayed_keys = set()

            for key in preferred_order:
                 if key in statistics and key != "Divergence Data": # Exclude divergence data from stats text area
                      value = statistics[key]
                      if isinstance(value, dict):
                          self.stats_text.insert(tk.END, f"{key}:\n")
                          for sub_key in sorted(value.keys()):
                               sub_value = value[sub_key]
                               if isinstance(sub_value, (float, int)):
                                    self.stats_text.insert(tk.END, f"  {sub_key}: {sub_value:.2f}\n" if isinstance(sub_value, float) else f"  {sub_key}: {sub_value}\n")
                               else:
                                     self.stats_text.insert(tk.END, f"  {sub_key}: {sub_value}\n")
                      elif isinstance(value, float):
                           self.stats_text.insert(tk.END, f"{key}: {value:.2f}\n")
                      elif isinstance(value, int):
                            self.stats_text.insert(tk.END, f"{key}: {value}\n")
                      else:
                           self.stats_text.insert(tk.END, f"{key}: {value}\n")
                      displayed_keys.add(key)

            # Add any other stats not in preferred_order
            for key, value in statistics.items():
                if key not in displayed_keys and key != "Divergence Data":
                     if isinstance(value, dict):
                          self.stats_text.insert(tk.END, f"{key}:\n")
                          for sub_key in sorted(value.keys()):
                               sub_value = value[sub_key]
                               if isinstance(sub_value, (float, int)):
                                    self.stats_text.insert(tk.END, f"  {sub_key}: {sub_value:.2f}\n" if isinstance(sub_value, float) else f"  {sub_key}: {sub_value}\n")
                               else:
                                     self.stats_text.insert(tk.END, f"  {sub_key}: {sub_value}\n")
                     elif isinstance(value, float):
                          self.stats_text.insert(tk.END, f"{key}: {value:.2f}\n")
                     elif isinstance(value, int):
                           self.stats_text.insert(tk.END, f"{key}: {value}\n")
                     else:
                          self.stats_text.insert(tk.END, f"{key}: {value}\n")
        else:
            self.stats_text.insert(tk.END, "No statistics available from backend.\n")
        self.stats_text.config(state="disabled")


        # --- Display Partition (NEXUS format example) ---
        self.partition_text.config(state="normal")
        self.partition_text.delete("1.0", tk.END)
        self.partition_text.insert(tk.END, "Partition (NEXUS format example):\n\n")
        if partition and self._concatenated_sequences: # Need sequences to get NTAX/NCHAR
             self.partition_text.insert(tk.END, "#nexus\n\n")

             num_taxa = len(self._concatenated_sequences)
             seq_len = len(next(iter(self._concatenated_sequences.values()))) if self._concatenated_sequences else 0

             if num_taxa > 0 and seq_len > 0:
                  self.partition_text.insert(tk.END, "BEGIN DATA;\n")
                  self.partition_text.insert(tk.END, f"  DIMENSIONS NTAX={num_taxa} NCHAR={seq_len};\n")
                  datatypes = set(item[2] for item in partition if item[2] != 'Unknown')
                  if 'DNA' in datatypes and 'Protein' in datatypes:
                       # Corrected f.write to use self.partition_text.insert
                       self.partition_text.insert(tk.END, "  FORMAT DATATYPE=DNA MISSING=- GAP=-;\n  [ WARNING: Mixed DNA/Protein data, DATATYPE=DNA may be unsuitable. ]\n")
                  elif 'Protein' in datatypes:
                       self.partition_text.insert(tk.END, "  FORMAT DATATYPE=Protein MISSING=- GAP=-;\n")
                  else:
                       self.partition_text.insert(tk.END, "  FORMAT DATATYPE=DNA MISSING=- GAP=-;\n  [ DATATYPE assumed as DNA ]\n")
                  self.partition_text.insert(tk.END, "END; [DATA]\n\n")
             else:
                  self.partition_text.insert(tk.END, "[ Cannot generate DATA block header - no concatenated sequences or zero length ]\n\n")

             self.partition_text.insert(tk.END, "BEGIN PAUP;\n")
             self.partition_text.insert(tk.END, "  [ Charsets for gene partitions ]\n")
             for backend_gene_name, coord_range, gene_type in partition:
                 clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name} = {coord_range};\n")

             self.partition_text.insert(tk.END, "\n  [ Charsets for codon positions (DNA genes only) ]\n")
             dna_genes_with_partitions = []
             for backend_gene_name, coord_range, gene_type in partition:
                  if gene_type and gene_type.upper() == 'DNA':
                       clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                       try:
                            start, end = map(int, coord_range.split('-'))
                            length = end - start + 1
                            if length >= 3:
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n")
                                 dna_genes_with_partitions.append(clean_gene_name)
                            elif length > 0:
                                 self.partition_text.insert(tk.END, f"  [ Warning: Gene {clean_gene_name} ({coord_range}) is too short ({length} bp) for full codon partitions. ]\n")
                       except ValueError:
                            self.partition_text.insert(tk.END, f"  [ Warning: Could not parse range '{coord_range}' for codon partitions for gene {clean_gene_name}. ]\n")
                  elif gene_type and gene_type.upper() == 'PROTEIN':
                       self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")
                  elif gene_type == 'Unknown':
                       self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} has Unknown type, no codon partitions generated. ]\n")


             if not dna_genes_with_partitions:
                  self.partition_text.insert(tk.END, "  [ No DNA genes found or genes too short for codon partitions. ]\n")


             self.partition_text.insert(tk.END, "\n  [ Link block (gene level) ]\n")
             if partition:
                 self.partition_text.insert(tk.END, "  link characters = ")
                 link_strings = []
                 for backend_gene_name, coord_range, gene_type in partition:
                      clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                      link_strings.append(f"{clean_gene_name} : {coord_range}")
                 self.partition_text.insert(tk.END, ", ".join(link_strings) + ";\n")
             else:
                  self.partition_text.insert(tk.END, "  [ No genes to link ]\n")


             # MrBayes Partition block
             self.partition_text.insert(tk.END, "\nBEGIN mrbayes;\n")
             self.partition_text.insert(tk.END, "  [ Example MrBayes partition block ]\n")
             self.partition_text.insert(tk.END, "  [ Gene level charsets ]\n")
             gene_charset_names_mb = []
             for backend_gene_name, coord_range, gene_type in partition:
                 clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name} = {coord_range};\n")
                 gene_charset_names_mb.append(clean_gene_name)

             self.partition_text.insert(tk.END, "\n  [ Codon position charsets (DNA genes only) ]\n")
             mr_bayes_dna_charsets_exist_mb = False
             pos1_names_mb = []
             pos2_names_mb = []
             pos3_names_mb = []

             for backend_gene_name, coord_range, gene_type in partition:
                  if gene_type and gene_type.upper() == 'DNA':
                       clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                       try:
                            start, end = map(int, coord_range.split('-'))
                            length = end - start + 1
                            if length >= 3:
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                 self.partition_text.insert(tk.END, f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n")
                                 mr_bayes_dna_charsets_exist_mb = True
                                 pos1_names_mb.append(f"{clean_gene_name}_pos1")
                                 pos2_names_mb.append(f"{clean_gene_name}_pos2")
                                 pos3_names_mb.append(f"{clean_gene_name}_pos3")
                       except ValueError:
                               pass

                  elif gene_type and gene_type.upper() == 'PROTEIN':
                        self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")
                  elif gene_type == 'Unknown':
                       self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} has Unknown type, no codon partitions generated. ]\n")


             if not mr_bayes_dna_charsets_exist_mb:
                  self.partition_text.insert(tk.END, "  [ No DNA genes found or genes too short for codon partitions. ]\n")

             self.partition_text.insert(tk.END, "\n  [ MrBayes Partition Commands (Choose ONE) ]\n")
             if gene_charset_names_mb:
                 self.partition_text.insert(tk.END, "  partition by_gene = {0}: {1};\n".format(
                     len(gene_charset_names_mb), " ".join(gene_charset_names_mb))
                 )
             if pos1_names_mb or pos2_names_mb or pos3_names_mb:
                  partition_groups_mrbayes_mb = []
                  if pos1_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos1_names_mb))
                  if pos2_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos2_names_mb))
                  if pos3_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos3_names_mb))
                  if partition_groups_mrbayes_mb:
                       total_mrbayes_partitions = len(partition_groups_mrbayes_mb)
                       self.partition_text.insert(tk.END, "  partition by_codon_pos = {0}: {1};\n".format(
                               total_mrbayes_partitions, ", ".join(partition_groups_mrbayes_mb))
                       )
                  else:
                       self.partition_text.insert(tk.END, "  [ Cannot generate 'by_codon_pos' partition command - check DNA charsets ]\n")

             self.partition_text.insert(tk.END, "\n  [ Example: To use 'by_gene' partition: ]\n")
             if gene_charset_names_mb:
                 self.partition_text.insert(tk.END, "  set partition = by_gene;\n")
             elif pos1_names_mb or pos2_names_mb or pos3_names_mb:
                  self.partition_text.insert(tk.END, "  [ Example: To use 'by_codon_pos' partition: ]\n")
                  self.partition_text.insert(tk.END, "  set partition = by_codon_pos;\n")
             else:
                 self.partition_text.insert(tk.END, "  [ No MrBayes partitions generated ]\n")

             self.partition_text.insert(tk.END, "END; [mrbayes]\n")

        else:
             self.partition_text.insert(tk.END, "No partition data available from backend or no concatenated sequences.\n")
        self.partition_text.config(state="disabled")


    # New method to update ONLY the divergence table display
    def _update_divergence_table_display(self, divergence_data, reference_taxon):
        """Updates the divergence table Treeview with data."""

        # Clear existing columns and data
        if self.divergence_treeview["columns"]:
             self.divergence_treeview["columns"] = ()

        self.divergence_treeview.delete(*self.divergence_treeview.get_children())
        self._treeview_id_to_taxon_name = {} # Clear map for new results

        if divergence_data:
             # Get columns from the keys of the first taxon's data
             sample_taxon_data = next(iter(divergence_data.values()), {})
             # Ensure 'Taxon', 'Total score' and 'No of charsets' are always first if present
             dynamic_columns = ["Taxon"]
             if "Total score" in sample_taxon_data: dynamic_columns.append("Total score")
             if "No of charsets" in sample_taxon_data: dynamic_columns.append("No of charsets")

             # Add gene columns, sorted alphabetically, excluding the standard keys
             gene_columns = sorted([key for key in sample_taxon_data.keys() if key not in ["Total score", "No of charsets"]])
             dynamic_columns.extend(gene_columns)


             self.divergence_treeview["columns"] = dynamic_columns
             for col in dynamic_columns:
                  heading_text = col
                  column_width = 100 # Default width
                  min_width = 50 # Default minwidth (safer than 0)
                  anchor = 'center'
                  stretch = tk.NO

                  if col == "Taxon":
                       heading_text = col
                       # Increased initial width and minwidth for Taxon column
                       column_width = 200
                       min_width = 150
                       anchor = 'w'
                       stretch = tk.YES # Ensure this column stretches
                  elif col in ["Total score", "No of charsets"]:
                       heading_text = col
                       column_width = 90
                       min_width = 80 # Set minwidth close to width
                       anchor = 'e'
                       stretch = tk.NO
                  else: # Assume it's a gene column
                        heading_text = col
                        column_width = 120
                        min_width = 100 # Set minwidth
                        anchor = 'center'
                        stretch = tk.NO

                  self.divergence_treeview.heading(col, text=heading_text, anchor=anchor)
                  # Apply width, minwidth, and stretch
                  self.divergence_treeview.column(col, width=column_width, minwidth=min_width, anchor=anchor, stretch=stretch)


             # Populate the table rows
             # Sort taxons alphabetically, but put the reference taxon first
             sorted_taxons = sorted(divergence_data.keys())
             if reference_taxon in sorted_taxons:
                 # Remove the reference taxon and insert it at the beginning
                 sorted_taxons.remove(reference_taxon)
                 sorted_taxons.insert(0, reference_taxon)


             reference_item_id = None # Variable to store the item ID of the reference taxon

             for taxon in sorted_taxons:
                  taxon_data = divergence_data.get(taxon, {}) # Use .get for safety
                  values = []

                  for col in dynamic_columns:
                       if col == "Taxon":
                            values.append(taxon)
                       elif col == "Total score":
                            score = taxon_data.get(col, 0.0)
                            if isinstance(score, (int, float)):
                                values.append(f"{score:.2f}%")
                            else:
                                values.append(str(score))

                       elif col == "No of charsets":
                             # Ensure default is 0 if key missing
                            values.append(taxon_data.get(col, 0))
                       else: # Gene columns
                            values.append(taxon_data.get(col, "N/A"))


                  # Use the current taxon name as the unique item ID.
                  # When the table is rebuilt, the item IDs will be the *current* names.
                  item_id = taxon

                  tags = ()
                  if taxon == reference_taxon:
                       tags = ('outgroup',) # Add 'outgroup' tag

                  # Insert the row with the determined item_id, values, and tags
                  inserted_id = self.divergence_treeview.insert("", tk.END, iid=item_id, values=values, tags=tags)
                  # Store the mapping from the Treeview item ID (which is the current taxon name string) to the current taxon name.
                  # This map is mainly used by the _on_treeview_click and _update_taxon_name logic.
                  self._treeview_id_to_taxon_name[inserted_id] = taxon

                  # If this is the reference taxon, store its new item ID (which is just the taxon name)
                  if taxon == reference_taxon:
                       reference_item_id = inserted_id


             # After populating, select the reference taxon row if it exists
             if reference_item_id:
                  # Clear any prior selection before setting the new one
                  current_selection = self.divergence_treeview.selection()
                  if current_selection:
                      self.divergence_treeview.selection_remove(*current_selection)
                  self.divergence_treeview.selection_set(reference_item_id)
                  # Optionally scroll to the selected item if the table is large
                  # self.divergence_treeview.see(reference_item_id)


        # No data case: Ensure headers are cleared and columns removed
        # Check if columns exist *after* attempting to set them
        if not self.divergence_treeview["columns"]:
             # Explicitly clear the main heading column (#0) if it's somehow visible without other columns
             try:
                  self.divergence_treeview.heading("#0", text="")
                  self.divergence_treeview.column("#0", width=0, minwidth=0, stretch=tk.NO)
             except tk.TclError:
                  # This might fail if #0 wasn't there in the first place, safe to ignore
                  pass

             # Ensure all potential numbered columns are also cleared
             # This loop might not be necessary if self.divergence_treeview["columns"] is already empty
             # but it adds safety if columns were set but then data was empty.
             for col_id in self.divergence_treeview.get_children(''): # This gets top-level item IDs, not column IDs. Incorrect usage.
                  # Correct way to get column identifiers if needed, but we just rely on `self.divergence_treeview["columns"] = ()`
                  pass # No action needed here

             # The line below is the primary way to clear defined columns
             self.divergence_treeview["columns"] = ()


    # New method to trigger recalculation and display of divergence
    def _recalculate_and_display_divergence(self):
        """
        Recalculates divergence based on the current _reference_taxon and updates the display.
        This uses the stored SequenceConcatenator instance and its internal data.
        """
        # Ensure the backend instance and necessary data exist
        if not self._concatenator or not self._concatenated_sequences or not self._backend_gene_info:
            # Should not happen if on_submit ran successfully and data was stored
            messagebox.showwarning("Recalculation Failed", "Cannot recalculate divergence. Necessary data is missing. Please re-run 'Concatenate & Analyze'.")
            return

        # Ensure the current reference taxon is still valid in the current concatenated data
        taxons_in_concat = sorted(list(self._concatenated_sequences.keys()))
        if not taxons_in_concat:
             # If no taxons are left, clear everything and show warning
             self._reference_taxon = None
             self._divergence_data = {}
             self._update_divergence_table_display(self._divergence_data, self._reference_taxon)
             messagebox.showwarning("No Taxons Available", "No taxons left to calculate divergence.")
             return

        if self._reference_taxon not in taxons_in_concat:
             # If the old reference taxon was removed or doesn't exist, revert to the first taxon
             old_ref = self._reference_taxon
             self._reference_taxon = taxons_in_concat[0] # Set to the new first taxon
             if old_ref: # Only show message if there was a previous reference
                 messagebox.showinfo("Reference Taxon Changed", f"Previous reference taxon '{old_ref}' is no longer available or not in current data. Setting reference to the first taxon: '{self._reference_taxon}'.")
             # If old_ref was None, it means there were taxons but no initial ref (unlikely after on_submit fix)


        try:
             # Call the public method on the stored backend instance to recalculate divergence
             # This method uses the data *within* the backend instance, only needing the reference taxon name
             # The backend handles updating its internal divergence data based on the new ref taxon.
             self._divergence_data = self._concatenator.recalculate_divergence_using_internal_data(self._reference_taxon)

             # Update the divergence table display with the new data and the current reference taxon
             self._update_divergence_table_display(self._divergence_data, self._reference_taxon)

        except Exception as e:
            messagebox.showerror("Divergence Recalculation Error", f"An error occurred while recalculating divergence:\n{e}")
            traceback.print_exc()
            # Optionally clear the divergence table on error
            self._divergence_data = {}
            self._update_divergence_table_display(self._divergence_data, self._reference_taxon)


    # New method triggered by the "Make Outgroup" button
    def make_selected_taxon_outgroup(self):
         """Sets the selected taxon in the divergence table as the outgroup (reference taxon)."""
         selected_items = self.divergence_treeview.selection()
         if not selected_items:
              messagebox.showwarning("Warning", "Please select a taxon in the table to make it the outgroup.")
              return
         if len(selected_items) > 1:
              messagebox.showwarning("Warning", "Please select only one taxon to make it the outgroup.")
              return

         selected_item_id = selected_items[0]
         # Get the taxon name using the map, which stores the current name associated with the item ID
         # With iid being the current name, the map should just return the iid itself.
         selected_taxon_name = self._treeview_id_to_taxon_name.get(selected_item_id, None)

         # Add a fallback to get the name directly from the Treeview if the map fails for some reason
         if selected_taxon_name is None:
              try:
                   # Assuming the taxon name is always the first value in the row
                   item_values = self.divergence_treeview.item(selected_item_id, 'values')
                   if item_values:
                        selected_taxon_name = item_values[0]
                        print(f"Warning: _treeview_id_to_taxon_name map missing entry for item '{selected_item_id}'. Using Treeview value '{selected_taxon_name}' as fallback.", file=sys.stderr)
                   else:
                        raise IndexError("No values found for item.")
              except (IndexError, tk.TclError):
                   messagebox.showerror("Internal Error", "Could not retrieve selected taxon name from selected item.")
                   return


         # Ensure the selected taxon name is actually present in the concatenated data keys
         if not self._concatenated_sequences or selected_taxon_name not in self._concatenated_sequences:
              messagebox.showwarning("Invalid Selection", f"Selected taxon '{selected_taxon_name}' is not a valid taxon in the concatenated data.")
              # Clear selection as it's invalid
              self.divergence_treeview.selection_remove(selected_item_id)
              return

         # Check if it's already the outgroup
         if selected_taxon_name == self._reference_taxon:
             # Clear selection but no recalculation needed
             self.divergence_treeview.selection_remove(selected_item_id)
             # messagebox.showinfo("Outgroup Selected", f"'{selected_taxon_name}' is already the current outgroup.") # Optional: inform user
             return

         # Set the new reference taxon
         self._reference_taxon = selected_taxon_name

         # Trigger recalculation and display update
         self._recalculate_and_display_divergence()


    # Kept as it allows editing taxon names directly in the results table, updated to handle outgroup tag
    def _on_treeview_click(self, event):
        # If there's an active editor, try to apply its change before starting a new one
        if self.entry_editor:
            self._update_taxon_name()

        # Identify the item and column clicked
        item_id = self.divergence_treeview.identify_row(event.y)
        col_identifier = self.divergence_treeview.identify_column(event.x)

        if not item_id or not col_identifier:
            return # Click was not on a valid item or column


        # Find the column index from its identifier (e.g., '#1', '#2')
        try:
             # Treeview columns start at #1 for the first data column
             col_index_in_columns_tuple = int(col_identifier.replace('#', '')) - 1
         # Catch ValueError if it's the implicit #0 column or some other identifier
        except ValueError:
             return # Not a standard data column identifier, ignore click for editing


        # Check if the first column ('Taxon' column) was clicked
        # Ensure columns exist first
        if self.divergence_treeview["columns"] and col_index_in_columns_tuple >= 0 and col_index_in_columns_tuple < len(self.divergence_treeview["columns"]):
             clicked_column_name = self.divergence_treeview["columns"][col_index_in_columns_tuple]

             if clicked_column_name == "Taxon":
                  bbox = self.divergence_treeview.bbox(item_id, col_identifier)
                  if not bbox: return

                  # Get the current taxon name associated with this item ID from the map
                  current_value = self._treeview_id_to_taxon_name.get(item_id, None)

                  # Fallback if map lookup fails (shouldn't happen with current logic, but safer)
                  if current_value is None:
                       try:
                            current_value = self.divergence_treeview.item(item_id, 'values')[col_index_in_columns_tuple]
                            print(f"Warning: _treeview_id_to_taxon_name map lookup failed for iid '{item_id}'. Using Treeview value '{current_value}' as fallback.", file=sys.stderr)
                       except (IndexError, tk.TclError):
                            print(f"Error: Could not get value from Treeview item '{item_id}', column '{col_identifier}'. Cannot start editor.", file=sys.stderr)
                            return # Cannot proceed

                  self.entry_item = item_id # Store the item ID
                  self.entry_column = col_index_in_columns_tuple # Store the column index

                  self._create_entry_widget(bbox, current_value)


    # Kept as it's part of the taxon renaming functionality
    def _create_entry_widget(self, bbox, initial_value):
        x, y, width, height = bbox

        # Destroy any existing editor before creating a new one
        if self.entry_editor:
             self._update_taxon_name() # Attempt to save changes of the previous editor

        self.entry_editor = tk.Entry(self.divergence_treeview, font=self.font_normal)
        self.entry_editor.insert(0, initial_value)

        # Bind events
        self.entry_editor.bind("<Return>", lambda event: self._update_taxon_name())
        self.entry_editor.bind("<Escape>", lambda event: self._cancel_edit())
        # Use a short delay for FocusOut to allow clicking buttons/other widgets before losing focus
        # However, binding FocusOut can sometimes interfere with other clicks.
        # A simple FocusOut bind is usually okay if it just saves. Let's keep the direct bind.
        self.entry_editor.bind("<FocusOut>", lambda event: self._update_taxon_name())


        # Place the editor widget
        # Need to get the absolute coordinates of the Treeview widget to place the Entry widget correctly on top of it
        treeview_x = self.divergence_treeview.winfo_rootx()
        treeview_y = self.divergence_treeview.winfo_rooty()
        widget_x = treeview_x + x
        widget_y = treeview_y + y

        # Place relative to the Treeview widget itself
        self.entry_editor.place(in_=self.divergence_treeview, x=x, y=y, width=width, height=height)

        self.entry_editor.focus_set()
        self.entry_editor.select_range(0, tk.END)


    # Kept as it's part of the taxon renaming functionality
    def _cancel_edit(self):
         """Destroys the entry editor without saving."""
         if self.entry_editor:
              self.entry_editor.destroy()
              self.entry_editor = None
              self.entry_item = None
              self.entry_column = None

    # Kept as it's part of the taxon renaming functionality, updated to handle outgroup tag
    def _update_taxon_name(self):
        if not self.entry_editor or self.entry_item is None or self.entry_column is None:
            return # No active editor or missing state

        new_name = self.entry_editor.get().strip()
        item_id_being_edited = self.entry_item
        col_index = self.entry_column

        # Get the old name using the stored item_id from the map
        old_name = self._treeview_id_to_taxon_name.get(item_id_being_edited)

        # Destroy the editor immediately
        self.entry_editor.destroy()
        self.entry_editor = None
        self.entry_item = None # Clear state after processing
        self.entry_column = None # Clear state after processing


        if not old_name:
             print(f"Error: Could not retrieve old taxon name from map for item ID '{item_id_being_edited}'. Aborting rename.", file=sys.stderr)
             messagebox.showerror("Internal Error", "Could not identify the taxon being edited.")
             return

        # If the name is empty or unchanged, just exit after destroying the editor
        if not new_name:
            messagebox.showwarning("Invalid Name", "Taxon name cannot be empty.")
            return
        if new_name == old_name:
            return

        # Check for duplicate names against the current keys in concatenated sequences (the source of truth)
        current_taxon_names = set(self._concatenated_sequences.keys()) if self._concatenated_sequences else set()
        # Temporarily remove the old name from the set for the duplicate check
        if old_name in current_taxon_names:
             current_taxon_names.remove(old_name)

        if new_name in current_taxon_names:
            messagebox.showwarning("Duplicate Name", f"A taxon named '{new_name}' already exists.")
            return

        # --- Update Internal Data Structures ---
        # Update concatenated sequences dictionary keys
        if self._concatenated_sequences and old_name in self._concatenated_sequences:
             seq = self._concatenated_sequences.pop(old_name)
             self._concatenated_sequences[new_name] = seq
             # Update __all_taxons in the backend instance if it exists, as it derives from this
             # WARNING: Directly accessing internal backend attributes is risky.
             # A public method in SequenceConcatenator to signal a rename would be better.
             if self._concatenator:
                  try:
                       current_backend_taxons = sorted(list(self._concatenated_sequences.keys()))
                       self._concatenator._SequenceConcatenator__all_taxons = current_backend_taxons
                       # Also update the taxon order which might be cached internally
                       if hasattr(self._concatenator, '_taxon_order'):
                           self._concatenator._taxon_order = current_backend_taxons
                  except AttributeError:
                       print("Warning: Could not update internal backend taxon list/order.", file=sys.stderr)


        # Update divergence data dictionary keys
        # Divergence data should be recalculated, but the dictionary structure itself needs the new key
        if self._divergence_data and old_name in self._divergence_data:
             # Save the data associated with the old name
             data = self._divergence_data.pop(old_name)
             # Add the data back using the new name as the key
             self._divergence_data[new_name] = data


        # Also update the reference taxon name if the renamed taxon was the reference
        if self._reference_taxon == old_name:
             self._reference_taxon = new_name


        # --- Update Treeview and Map ---
        # The easiest way to update the Treeview display and the internal map
        # after keys in _concatenated_sequences and _divergence_data have changed
        # is to clear and repopulate the Treeview. This is handled by
        # _update_divergence_table_display, which rebuilds based on the
        # *current* state of _divergence_data and sets item IDs to the *current* names.
        # The map _treeview_id_to_taxon_name is cleared and rebuilt within this function.

        # Recalculate divergence if the *reference* taxon name changed.
        # If a non-reference taxon was renamed, the divergence values relative
        # to the current reference don't change, but the row needs to show the new name.
        # Calling _recalculate_and_display_divergence handles both cases:
        # 1. If _reference_taxon changed, it recalculates values via the backend.
        # 2. It always calls _update_divergence_table_display which rebuilds the UI
        #    and the map using the latest data (with the new name).
        self._recalculate_and_display_divergence()

        # No explicit treeview item update needed here because _recalculate_and_display_divergence
        # (which calls _update_divergence_table_display) rebuilds the table.


    # Removed as the toolbar buttons that called it are removed
    # def remove_selected_taxon_from_results(self): ...
    # Removed as the toolbar buttons that called it are removed
    # def excise_entire_taxon_from_results(self): ...


    # Kept as it's used internally by other kept functions (remove_selected_gene, on_submit, on_reset)
    def clear_results_display(self):
        """Clears the text areas and treeview in the Results tab and stored results."""

        for text_widget in [self.stats_text, self.partition_text]:
            text_widget.config(state="normal")
            text_widget.delete("1.0", tk.END)
            text_widget.config(state="disabled")

        # Clear Treeview columns and data
        if self.divergence_treeview["columns"]:
             self.divergence_treeview["columns"] = ()
        self.divergence_treeview.delete(*self.divergence_treeview.get_children())
        self._treeview_id_to_taxon_name = {} # Clear map

        # Clear stored results
        self._concatenated_sequences = None
        self._concatenator = None # Clear the backend instance
        self._backend_gene_info = None
        self._partition_data = None
        self._statistics = None
        self._divergence_data = None
        self._reference_taxon = None


    # Kept as it's triggered by the "Reset All" button on the Input tab
    def on_reset(self):
        """Resets the input (genes list) and clears results."""
        # Check if there is anything to reset before asking
        if self.gene_names or self.raw_gene_contents or self._concatenated_sequences is not None or self._divergence_data is not None:
            confirm = messagebox.askyesno("Confirm Reset", "Are you sure you want to remove all loaded genes and clear results?")
            if not confirm:
                return
        else:
             # Nothing to reset, just return
             return


        self.gene_names = []
        self.raw_gene_contents = []
        self.gene_listbox.delete(0, tk.END)

        self.clear_results_display()

    # Restored as requested
    def export_fasta(self):
        """Exports the concatenated sequences to a FASTA file using current taxon names."""
        if not self._concatenated_sequences:
            messagebox.showwarning("Warning", "No concatenated sequences to export. Run 'Concatenate & Analyze' first.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=(("FASTA files", "*.fasta *.fas *.fna *.faa *.fa"), ("All files", "*.*")),
            title="Save Concatenated FASTA File"
        )

        if not filepath:
            return

        try:
            with open(filepath, "w", encoding='utf-8') as f:
                # Use the keys from the potentially edited _concatenated_sequences
                sorted_taxons = sorted(self._concatenated_sequences.keys())
                for taxon in sorted_taxons:
                    sequence = self._concatenated_sequences[taxon]
                    # Simple cleaning for FASTA header (remove leading/trailing whitespace)
                    # Also replace any characters that might cause issues in FASTA header lines like >
                    clean_taxon_header = re.sub(r'[>\s]', '_', taxon.strip()) # Replace > and whitespace with underscore
                    f.write(f">{clean_taxon_header}\n")
                    wrapped_seq = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
                    f.write(f"{wrapped_seq}\n")

            messagebox.showinfo("Export Successful", f"Concatenated sequences exported to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save FASTA file:\n{e}")
            traceback.print_exc()

    # Restored as requested
    def export_partition(self):
        """Exports the partition data to a NEXUS formatted file (partition block only)."""
        # Need partition data AND concatenated sequences (for NTAX/NCHAR)
        if not self._partition_data or not self._concatenated_sequences:
            messagebox.showwarning("Warning", "No partition data or concatenated sequences available. Run 'Concatenate & Analyze' first.")
            return

        num_taxa = len(self._concatenated_sequences)
        seq_len = len(next(iter(self._concatenated_sequences.values()))) if self._concatenated_sequences else 0

        if num_taxa == 0 or seq_len == 0:
             messagebox.showwarning("Warning", "Concatenated sequence data is needed to determine matrix dimensions for the partition file header. Run 'Concatenate & Analyze' successfully first.")
             return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".nexus",
            filetypes=(("NEXUS files", "*.nexus *.nex"), ("Text files", "*.txt"), ("All files", "*.*")),
            title="Save Partition File (NEXUS partition block)"
        )

        if not filepath:
            return

        try:
            with open(filepath, "w", encoding='utf-8') as f:
                f.write("#nexus\n\n")

                f.write("BEGIN DATA;\n")
                f.write(f"  DIMENSIONS NTAX={num_taxa} NCHAR={seq_len};\n")
                datatypes = set(item[2] for item in self._partition_data if item[2] != 'Unknown')
                if 'DNA' in datatypes and 'Protein' in datatypes:
                     f.write("  FORMAT DATATYPE=DNA MISSING=- GAP=-;\n  [ WARNING: Mixed DNA/Protein data, DATATYPE=DNA may be unsuitable. ]\n")
                elif 'Protein' in datatypes:
                     f.write("  FORMAT DATATYPE=Protein MISSING=- GAP=-;\n")
                else:
                     f.write("  FORMAT DATATYPE=DNA MISSING=- GAP=-;\n  [ DATATYPE assumed as DNA ]\n")
                f.write("END; [DATA]\n\n")

                f.write("BEGIN PAUP;\n")
                f.write("  [ Charsets for gene partitions ]\n")
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                    clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                    f.write(f"  charset {clean_gene_name} = {coord_range};\n")

                f.write("\n  [ Charsets for codon positions (DNA genes only) ]\n")
                dna_genes_with_partitions_in_partition = []
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                     if gene_type and gene_type.upper() == 'DNA':
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
                     elif gene_type and gene_type.upper() == 'PROTEIN':
                          f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")
                     elif gene_type == 'Unknown':
                          f.write(f"  [ Note: Gene {backend_gene_name} has Unknown type, no codon partitions generated. ]\n")


                if not dna_genes_with_partitions_in_partition:
                     f.write("  [ No DNA genes found or genes too short for codon partitions. ]\n")

                f.write("\n  [ Link block (gene level) ]\n")
                if self._partition_data:
                     f.write("  link characters = ")
                     link_strings = []
                     for backend_gene_name, coord_range, gene_type in self._partition_data:
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          link_strings.append(f"{clean_gene_name} : {coord_range}")
                     f.write(", ".join(link_strings) + ";\n")
                else:
                     f.write("  [ No genes to link ]\n")

                # MrBayes Partition block for export
                f.write("\nBEGIN mrbayes;\n")
                f.write("  [ Example MrBayes partition block ]\n")
                f.write("  [ Gene level charsets ]\n")
                gene_charset_names_mb = []
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
                     if gene_type and gene_type.upper() == 'DNA':
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
                          except ValueError:
                               pass

                     elif gene_type and gene_type.upper() == 'PROTEIN':
                        f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")
                     elif gene_type == 'Unknown':
                          f.write(f"  [ Note: Gene {backend_gene_name} has Unknown type, no codon partitions generated. ]\n")


                if not mr_bayes_dna_charsets_exist_mb:
                     f.write("  [ No DNA genes found or genes too short for codon partitions. ]\n")

                f.write("\n  [ MrBayes Partition Commands (Choose ONE) ]\n")
                if gene_charset_names_mb:
                    f.write("  partition by_gene = {0}: {1};\n".format(
                        len(gene_charset_names_mb), " ".join(gene_charset_names_mb))
                    )
                if pos1_names_mb or pos2_names_mb or pos3_names_mb:
                     partition_groups_mrbayes_mb = []
                     if pos1_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos1_names_mb))
                     if pos2_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos2_names_mb))
                     if pos3_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos3_names_mb))
                     if partition_groups_mrbayes_mb:
                          total_mrbayes_partitions = len(partition_groups_mrbayes_mb)
                          f.write("  partition by_codon_pos = {0}: {1};\n".format(
                               total_mrbayes_partitions, ", ".join(partition_groups_mrbayes_mb))
                          )
                     else:
                          f.write("  [ Cannot generate 'by_codon_pos' partition command - check DNA charsets ]\n")

                f.write("\n  [ Example: To use 'by_gene' partition: ]\n")
                if gene_charset_names_mb:
                    f.write("  set partition = by_gene;\n")
                elif pos1_names_mb or pos2_names_mb or pos3_names_mb:
                     f.write("  [ Example: To use 'by_codon_pos' partition: ]\n")
                     f.write("  set partition = by_codon_pos;\n")
                else:
                    f.write("  [ No MrBayes partitions generated ]\n")

                f.write("END; [mrbayes]\n")


            messagebox.showinfo("Export Successful", f"Partition data exported to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save partition file:\n{e}")
            traceback.print_exc()

    # Restored as requested
    def export_full_nexus(self):
        """Exports the concatenated sequences and partition data to a full NEXUS file using current taxon names."""
        if not self._concatenated_sequences or not self._partition_data:
            messagebox.showwarning("Warning", "No concatenated sequences or partition data available. Run 'Concatenate & Analyze' first.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".nexus",
            filetypes=(("NEXUS files", "*.nexus *.nex"), ("All files", "*.*")),
            title="Save Full NEXUS File"
        )

        if not filepath:
            return

        try:
            with open(filepath, "w", encoding='utf-8') as f:
                f.write("#nexus\n\n")

                sorted_taxons = sorted(self._concatenated_sequences.keys())
                num_taxa = len(sorted_taxons)
                seq_len = len(next(iter(self._concatenated_sequences.values()))) if self._concatenated_sequences else 0

                if num_taxa == 0 or seq_len == 0:
                     messagebox.showwarning("Warning", "No concatenated sequences or zero length. Cannot export NEXUS matrix.")
                     return # Do not write an incomplete file


                # TAXA block
                f.write("BEGIN TAXA;\n")
                f.write(f"  DIMENSIONS NTAX={num_taxa};\n")
                f.write("  TAXLABELS\n")

                # Quote taxon labels if they contain special characters or whitespace
                for taxon in sorted_taxons:
                    # Using a slightly more comprehensive regex for chars requiring quoting in NEXUS
                    if re.search(r'[\s\'"`=;:,\[\]\(\)]+', str(taxon)) or taxon.strip() == '':
                         f.write(f"    '{taxon}'\n")
                    else:
                         f.write(f"    {taxon}\n")
                f.write("  ;\n")
                f.write("END; [TAXA]\n\n")


                # DATA block (Matrix and Format)
                f.write("BEGIN DATA;\n")
                f.write(f"  DIMENSIONS NTAX={num_taxa} NCHAR={seq_len};\n")

                datatypes = set(item[2] for item in self._partition_data if item[2] != 'Unknown')
                if 'DNA' in datatypes and 'Protein' in datatypes:
                     f.write("  FORMAT DATATYPE=DNA MISSING=- GAP=- INTERLEAVE=NO;\n  [ WARNING: Mixed DNA/Protein data, DATATYPE=DNA may be unsuitable. ]\n")
                elif 'Protein' in datatypes:
                     f.write("  FORMAT DATATYPE=Protein MISSING=- GAP=- INTERLEAVE=NO;\n")
                else:
                     f.write("  FORMAT DATATYPE=DNA MISSING=- GAP=- INTERLEAVE=NO;\n  [ DATATYPE assumed as DNA ]\n")

                f.write("  MATRIX;\n")

                for taxon in sorted_taxons:
                    sequence = self._concatenated_sequences[taxon]
                    # Quote taxon labels as done in TAXA block
                    taxon_label = f"'{taxon}'" if re.search(r'[\s\'"`=;:,\[\]\(\)]+', str(taxon)) or taxon.strip() == '' else taxon

                    f.write(f"{taxon_label}\n")
                    wrapped_seq = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
                    # Indent wrapped sequence lines appropriately for NEXUS MATRIX format
                    indented_wrapped_seq = '\n'.join([f"  {line}" for line in wrapped_seq.splitlines()])
                    f.write(f"{indented_wrapped_seq}\n")


                f.write("  ;\n")
                f.write("END; [DATA]\n\n")

                # PAUP block (Charsets and Link)
                f.write("BEGIN PAUP;\n")
                f.write("  [ Charsets for gene partitions ]\n")
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                    clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                    f.write(f"  charset {clean_gene_name} = {coord_range};\n")

                f.write("\n  [ Charsets for codon positions (DNA genes only) ]\n")
                dna_genes_with_partitions_in_nexus = []
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                     if gene_type and gene_type.upper() == 'DNA':
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
                     elif gene_type and gene_type.upper() == 'PROTEIN':
                           f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")
                     elif gene_type == 'Unknown':
                          f.write(f"  [ Note: Gene {backend_gene_name} has Unknown type, no codon partitions generated. ]\n")


                if not dna_genes_with_partitions_in_nexus:
                     f.write("  [ No DNA genes found or genes too short for codon partitions. ]\n")

                f.write("\n  [ Link block (gene level) ]\n")
                if self._partition_data:
                    f.write("  link characters = ")
                    link_strings = []
                    for backend_gene_name, coord_range, gene_type in self._partition_data:
                         clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                         link_strings.append(f"{clean_gene_name} : {coord_range}")
                    f.write(", ".join(link_strings) + ";\n")
                else:
                    f.write("  [ No genes to link ]\n")

                f.write("END; [PAUP]\n\n")

                # MrBayes Partition block
                f.write("BEGIN mrbayes;\n")
                f.write("  [ Example MrBayes partition block ]\n")
                f.write("  [ Gene level charsets ]\n")
                gene_charset_names_mb = []
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
                     if gene_type and gene_type.upper() == 'DNA':
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
                          except ValueError:
                               pass

                     elif gene_type and gene_type.upper() == 'PROTEIN':
                        f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")
                     elif gene_type == 'Unknown':
                          f.write(f"  [ Note: Gene {backend_gene_name} has Unknown type, no codon partitions generated. ]\n")


                if not mr_bayes_dna_charsets_exist_mb:
                     f.write("  [ No DNA genes found or genes too short for codon partitions. ]\n")

                f.write("\n  [ MrBayes Partition Commands (Choose ONE) ]\n")
                if gene_charset_names_mb:
                    f.write("  partition by_gene = {0}: {1};\n".format(
                        len(gene_charset_names_mb), " ".join(gene_charset_names_mb))
                    )
                if pos1_names_mb or pos2_names_mb or pos3_names_mb:
                     partition_groups_mrbayes_mb = []
                     if pos1_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos1_names_mb))
                     if pos2_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos2_names_mb))
                     if pos3_names_mb: partition_groups_mrbayes_mb.append(" ".join(pos3_names_mb))
                     if partition_groups_mrbayes_mb:
                          total_mrbayes_partitions = len(partition_groups_mrbayes_mb)
                          f.write("  partition by_codon_pos = {0}: {1};\n".format(
                               total_mrbayes_partitions, ", ".join(partition_groups_mrbayes_mb))
                          )
                     else:
                          f.write("  [ Cannot generate 'by_codon_pos' partition command - check DNA charsets ]\n")

                f.write("\n  [ Example: To use 'by_gene' partition: ]\n")
                if gene_charset_names_mb:
                    f.write("  set partition = by_gene;\n")
                elif pos1_names_mb or pos2_names_mb or pos3_names_mb:
                     f.write("  [ Example: To use 'by_codon_pos' partition: ]\n")
                     f.write("  set partition = by_codon_pos;\n")
                else:
                    f.write("  [ No MrBayes partitions generated ]\n")

                f.write("END; [mrbayes]\n")


            messagebox.showinfo("Export Successful", f"Full NEXUS file exported to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save Full NEXUS file:\n{e}")
            traceback.print_exc()


# --- Main execution ---
if __name__ == "__main__":
    root = tk.Tk()
    app = SequenceConcatenatorApp(root)
    root.mainloop()