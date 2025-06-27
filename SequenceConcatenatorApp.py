import sys
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import re
import os
import traceback

# Assuming SequenceConcatenator.py is in the same directory or accessible in the python path
try:
    from SequenceConcatenator import SequenceConcatenator
except ImportError:
    messagebox.showerror("Import Error", "Could not find SequenceConcatenator.py. Please ensure it is in the same directory.")
    exit()


class SequenceConcatenatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Sequence Concatenator UI") # Simple title
        self.root.geometry("1200x900") # Increased size
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
        self.font_normal = ("Arial", 10) # Slightly smaller font for table
        self.font_bold = ("Arial", 10, "bold")
        self.font_heading = ("Arial", 11, "bold") # Slightly smaller heading
        self.font_button = ("Segoe UI", 10, "bold")
        self.font_mono = ("Courier New", 9) # Monospaced for sequence/partition display


        # --- Menu Bar (Minimal) ---
        menubar = tk.Menu(root)
        root.config(menu=menubar)

        filemenu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Exit", command=root.quit)

        importmenu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Import", menu=importmenu)
        importmenu.add_command(label="Import Genes...", command=self.open_add_gene_dialog)

        exportmenu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Export", menu=exportmenu)
        exportmenu.add_command(label="Export Concatenated FASTA...", command=self.export_fasta)
        exportmenu.add_command(label="Export Partition (NEXUS)...", command=self.export_partition)
        exportmenu.add_command(label="Export Full NEXUS...", command=self.export_full_nexus)


        # --- Toolbar (Minimal, only results-related actions requested) ---
        toolbar = tk.Frame(root, bg="#dddddd")
        toolbar.pack(side="top", fill="x")

        # Keep only taxon removal as requested implicitly by the image/description
        ttk.Button(toolbar, text="Delete taxon", command=lambda: self.remove_selected_taxon_from_results()).pack(side="left", padx=2, pady=2)
        ttk.Button(toolbar, text="Excise entire taxon", command=lambda: self.excise_entire_taxon_from_results()).pack(side="left", padx=2, pady=2)


        # Notebook (Tabs)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(pady=10, padx=10, fill="both", expand=True)

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
        style.map("Treeview", background=[("selected", self.primary_color)], foreground=[("selected", "white")])


        for tab_name, tab_frame in self.tabs.items():
            self.notebook.add(tab_frame, text=tab_name)

        self.setup_input_tab()

        # Container frame for results tab content
        self.results_content_frame = ttk.Frame(self.tabs["Results"], style='TFrame')
        self.results_content_frame.pack(fill="both", expand=True)
        self.setup_results_tab(self.results_content_frame)


        # Store results internally for export and display
        self._concatenated_sequences = None # {taxon_name: sequence_string, ...} -> Keys are potentially edited
        self._partition_data = None # [(gene_name, 'start-end', gene_type), ...] -> Uses original gene names from backend
        self._statistics = None # Dictionary holding various stats, including divergence
        self._divergence_data = None # {taxon_name: {'Total score': float, 'No of charsets': int, 'gene1_name': 'count (perc%)', ...}, ...}

        # Variable to hold the Treeview Entry widget for editing
        self.entry_editor = None
        self.entry_item = None
        self.entry_column = None
        # Map to track Treeview item IDs to current taxon names for editing/deleting
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
        partition_scrollbar_x.grid(row=1, column=0, sticky="ew") # Place below the text area
        self.partition_text.config(xscrollcommand=partition_scrollbar_x.set)


        # --- Divergence Table (Replaces Sequences Display) ---
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
        self.divergence_treeview.config(yscrollcommand=tree_scrollbar_y.set) # Fixed typo here as well

        tree_scrollbar_x = tk.Scrollbar(divergence_frame, orient="horizontal", command=self.divergence_treeview.xview)
        tree_scrollbar_x.grid(row=1, column=0, sticky="ew")
        self.divergence_treeview.config(xscrollcommand=tree_scrollbar_x.set)


        # Bind click for editing
        self.divergence_treeview.bind('<ButtonRelease-1>', self._on_treeview_click)


        # --- Export Buttons ---
        export_btn_frame = tk.Frame(parent_frame, bg=self.bg_color)
        export_btn_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10), padx=10)

        export_fasta_btn = tk.Button(
            export_btn_frame, text="Export FASTA", command=self.export_fasta,
            bg=self.button_bg_primary, fg=self.button_fg,
            activebackground=self.active_bg_primary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_fasta_btn.pack(side="left", padx=(0, 10))

        export_partition_btn = tk.Button(
            export_btn_frame, text="Export Partition (NEXUS block)", command=self.export_partition,
            bg=self.button_bg_primary, fg=self.button_fg,
            activebackground=self.active_bg_primary, activeforeground=self.button_fg,
            font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_partition_btn.pack(side="left", padx=(0, 10))

        export_full_nexus_btn = tk.Button(
             export_btn_frame, text="Export Full NEXUS", command=self.export_full_nexus,
             bg=self.button_bg_primary, fg=self.button_fg,
             activebackground=self.active_bg_primary, activeforeground=self.button_fg,
             font=self.font_button, bd=0, relief="flat", padx=20, pady=10, cursor="hand2"
        )
        export_full_nexus_btn.pack(side="left", padx=(0, 10))


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
                 dialog.destroy()
             else:
                  messagebox.showerror("Load Failed", status_msg, parent=dialog)
        else:
             pass


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

    def remove_selected_gene(self):
        selected_indices = self.gene_listbox.curselection()
        if not selected_indices:
            messagebox.showwarning("Warning", "Please select one or more genes to remove.")
            return

        selected_gene_names = [self.gene_listbox.get(i) for i in selected_indices]
        confirm = messagebox.askyesno("Confirm Removal", f"Are you sure you want to remove the selected gene(s)?\n\n" + "\n".join(selected_gene_names))
        if not confirm:
            return

        for i in sorted(selected_indices, reverse=True):
            self.gene_listbox.delete(i)
            del self.gene_names[i]
            del self.raw_gene_contents[i]

        self.clear_results_display()

    def on_submit(self):
        if not self.raw_gene_contents:
            messagebox.showwarning("Warning", "Please add at least one gene before submitting.")
            return

        self.clear_results_display()
        self._treeview_id_to_taxon_name = {} # Clear map for new results

        try:
            concatenator = SequenceConcatenator(self.raw_gene_contents)

            self._concatenated_sequences = concatenator.get_concatenated_sequences()
            self._statistics = concatenator.get_statistics()
            self._partition_data = concatenator.get_partition()
            self._divergence_data = self._statistics.get("Divergence Data", {})

            if not self._concatenated_sequences:
                 messagebox.showwarning("Concatenation Warning", "No concatenated sequences were produced by the backend. This might happen if no valid data was parsed or there are no common taxa across genes.")
                 self.clear_results_display()
                 return
            if not self._divergence_data and self._concatenated_sequences and len(self._concatenated_sequences) > 1: # Warn only if >1 taxon but no divergence
                 messagebox.showwarning("Divergence Warning", "No divergence data was produced by the backend. This might indicate an issue with calculating scores.")


            self.display_results(self._concatenated_sequences, self._statistics, self._partition_data, self._divergence_data)

            self.notebook.select(self.tabs["Results"])

        except Exception as e:
            messagebox.showerror("Processing Error", f"An error occurred during sequence processing in the backend.\nError details: {e}")
            traceback.print_exc()
            self.clear_results_display()


    def display_results(self, concatenated_sequences, statistics, partition, divergence_data):
        """Populates the Results tab widgets with data, including the Treeview."""

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
                "Missing Data per Taxon (count)",
                "Missing Data per Gene (count)",
                "Gene Lengths (in concatenated alignment)",
            ]
            displayed_keys = set()

            for key in preferred_order:
                 if key in statistics and key != "Divergence Data":
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
        if partition:
             self.partition_text.insert(tk.END, "#nexus\n\n")

             num_taxa = len(concatenated_sequences) if concatenated_sequences else 0
             seq_len = len(list(concatenated_sequences.values())[0]) if concatenated_sequences and list(concatenated_sequences.values()) else 0
             if num_taxa > 0 and seq_len > 0:
                  self.partition_text.insert(tk.END, "BEGIN DATA;\n")
                  self.partition_text.insert(tk.END, f"  DIMENSIONS NTAX={num_taxa} NCHAR={seq_len};\n")
                  datatypes = set(item[2] for item in partition if item[2] != 'Unknown')
                  if 'DNA' in datatypes and 'Protein' in datatypes:
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
                  if gene_type.upper() == 'DNA':
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
                  elif gene_type.upper() == 'PROTEIN':
                       self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")


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


             if dna_genes_with_partitions:
                self.partition_text.insert(tk.END, "\n  [ Link block (by codon position, DNA genes only) ]\n")
                self.partition_text.insert(tk.END, "  link characters = ")
                pos1_links = []
                pos2_links = []
                pos3_links = []

                for name, range_str, gene_type in partition:
                     if gene_type.upper() == 'DNA':
                          clean_name = re.sub(r'\W+', '_', name)
                          try:
                               start, end = map(int, range_str.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    pos1_links.append(f"{clean_name}_pos1 : {start}-{end}\\3")
                                    pos2_links.append(f"{clean_name}_pos2 : {start+1}-{end}\\3")
                                    pos3_links.append(f"{clean_name}_pos3 : {start+2}-{end}\\3")
                          except ValueError:
                                   pass

                all_codon_links = pos1_links + pos2_links + pos3_links
                if all_codon_links:
                    self.partition_text.insert(tk.END, ", ".join(all_codon_links) + ";\n")
                else:
                     self.partition_text.insert(tk.END, "[ No DNA genes long enough for codon link block ]\n")


             self.partition_text.insert(tk.END, "END; [PAUP]\n\n")

             self.partition_text.insert(tk.END, "BEGIN mrbayes;\n")
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
                  if gene_type.upper() == 'DNA':
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
                                 pos2_names_mb.append(f"{clean_name}_pos2") # Fix: Use clean_name_pos2
                                 pos3_names_mb.append(f"{clean_name}_pos3") # Fix: Use clean_name_pos3
                       except ValueError:
                               pass
                  elif gene_type.upper() == 'PROTEIN':
                        self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")

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
                       self.partition_text.insert(tk.END, "  partition by_codon_pos = {0}: {1};\n".format(
                            len(partition_groups_mrbayes_mb), ", ".join(partition_groups_mrbayes_mb))
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
                 self.partition_text.insert(tk.END, "  [ No partitions generated ]\n")

             self.partition_text.insert(tk.END, "END; [mrbayes]\n")

        else:
             self.partition_text.insert(tk.END, "No partition data available from backend.\n")
        self.partition_text.config(state="disabled")


        # --- Display Divergence Table ---
        # Clear existing columns and data
        # Need to get columns BEFORE deleting them from the Treeview object
        current_cols = list(self.divergence_treeview["columns"])
        for col in current_cols:
             # Must remove headings first for some reason with delete(col)
             self.divergence_treeview.heading(col, text="")
             self.divergence_treeview.column(col, width=0, stretch=tk.NO) # Also hide/shrink column
             self.divergence_treeview.delete(col) # Delete old columns - this line is incorrect usage, remove it

        # Correct way to clear columns:
        if self.divergence_treeview["columns"]: # Check if there are any columns to delete
             self.divergence_treeview["columns"] = () # Assign empty tuple to clear all columns


        self.divergence_treeview.delete(*self.divergence_treeview.get_children()) # Delete old rows
        self._treeview_id_to_taxon_name = {} # Clear map before repopulating

        if divergence_data:
             sample_taxon_data = next(iter(divergence_data.values()), {})
             # Ensure 'Total score' and 'No of charsets' are always first if present
             dynamic_columns = ["Taxon"]
             if "Total score" in sample_taxon_data: dynamic_columns.append("Total score")
             if "No of charsets" in sample_taxon_data: dynamic_columns.append("No of charsets")

             # Add gene columns, sorted alphabetically
             gene_columns = sorted([key for key in sample_taxon_data.keys() if key not in ["Total score", "No of charsets"]])
             dynamic_columns.extend(gene_columns)


             self.divergence_treeview["columns"] = dynamic_columns
             for col in dynamic_columns:
                  heading_text = col
                  column_width = 100
                  anchor = 'center'

                  if col == "Taxon":
                       column_width = 150
                       anchor = 'w'
                       stretch = tk.YES # Make Taxon column stretchable
                  elif col in ["Total score", "No of charsets"]:
                       column_width = 90
                       anchor = 'e'
                       stretch = tk.NO
                  else: # Assume it's a gene column
                        column_width = 120
                        anchor = 'center'
                        stretch = tk.NO


                  self.divergence_treeview.heading(col, text=heading_text, anchor=anchor)
                  self.divergence_treeview.column(col, width=column_width, anchor=anchor, stretch=stretch) # Use calculated stretch


             sorted_taxons = sorted(divergence_data.keys())
             for taxon in sorted_taxons:
                  taxon_data = divergence_data[taxon]
                  values = []

                  # Populate values list according to dynamic_columns order
                  for col in dynamic_columns:
                       if col == "Taxon":
                            values.append(taxon)
                       elif col == "Total score":
                            values.append(f"{taxon_data.get(col, 0.0):.2f}%")
                       elif col == "No of charsets":
                            values.append(taxon_data.get(col, 0))
                       else: # Gene columns
                            values.append(taxon_data.get(col, "N/A"))


                  # Insert the row, use the initial taxon name as the item ID
                  item_id = taxon # Use the taxon name as the unique item ID
                  self.divergence_treeview.insert("", tk.END, iid=item_id, values=values)
                  # Store the mapping from the Treeview item ID to the current taxon name
                  self._treeview_id_to_taxon_name[item_id] = taxon # Initially, item ID == taxon name


        else:
             self.divergence_treeview.heading("#0", text="")
             self.divergence_treeview.column("#0", width=0, stretch=tk.NO)
             self.divergence_treeview["columns"] = ()


    def _on_treeview_click(self, event):
        if self.entry_editor:
            self._update_taxon_name()

        item = self.divergence_treeview.identify_row(event.y)
        col_identifier = self.divergence_treeview.identify_column(event.x)

        # Check if a valid row and the first column ('#1' corresponding to the first column in "columns" list) were clicked
        if item and col_identifier == '#1': # '#1' is the first column displayed, which is 'Taxon'
            col_index = 0 # The 'Taxon' column is always the first column in our 'columns' list

            bbox = self.divergence_treeview.bbox(item, col_identifier)
            if not bbox: return

            current_value = self.divergence_treeview.item(item, 'values')[col_index]

            self.entry_item = item # Treeview item ID (original taxon name)
            self.entry_column = col_index

            self._create_entry_widget(bbox, current_value)


    def _create_entry_widget(self, bbox, initial_value):
        x, y, width, height = bbox

        self.entry_editor = tk.Entry(self.divergence_treeview, font=self.font_normal)
        self.entry_editor.insert(0, initial_value)

        self.entry_editor.bind("<Return>", lambda event: self._update_taxon_name())
        self.entry_editor.bind("<FocusOut>", lambda event: self._update_taxon_name())

        self.entry_editor.place(x=x, y=y, width=width, height=height, anchor='nw')
        self.entry_editor.focus_set()
        self.entry_editor.select_range(0, tk.END)


    def _update_taxon_name(self):
        if not self.entry_editor:
            return

        new_name = self.entry_editor.get().strip()
        item = self.entry_item
        col_index = self.entry_column # Should be 0

        # Get the current name using the map associated with this item ID
        # The map is the authoritative source for the current name associated with this Treeview item ID
        old_name = self._treeview_id_to_taxon_name.get(item)
        if old_name is None: # Should not happen if map is populated correctly
             print(f"Error: Map entry not found for item ID '{item}'. Using Treeview value as fallback.", file=sys.stderr)
             old_name = self.divergence_treeview.item(item, 'values')[col_index]


        if not new_name:
            messagebox.showwarning("Invalid Name", "Taxon name cannot be empty.")
            self.entry_editor.destroy()
            self.entry_editor = None
            self.entry_item = None
            self.entry_column = None
            return
        if new_name == old_name:
            self.entry_editor.destroy()
            self.entry_editor = None
            self.entry_item = None
            self.entry_column = None
            return

        # Check for duplicate names across *all* current taxon names in the internal data
        # We need to check against the *currently active* names, which are the keys in the internal dicts.
        # Let's use the keys from _concatenated_sequences as the source of truth for current names.
        current_taxon_names = set(self._concatenated_sequences.keys()) if self._concatenated_sequences else set()

        # Remove the old name from the set for the check (since we're renaming it)
        if old_name in current_taxon_names:
            current_taxon_names.remove(old_name)

        if new_name in current_taxon_names:
            messagebox.showwarning("Duplicate Name", f"A taxon named '{new_name}' already exists.")
            self.entry_editor.destroy()
            self.entry_editor = None
            self.entry_item = None
            self.entry_column = None
            return

        # Update internal data structures using the old_name (from map) and new_name
        # Update the key in the concatenated sequences dictionary
        if self._concatenated_sequences and old_name in self._concatenated_sequences:
             seq = self._concatenated_sequences.pop(old_name)
             self._concatenated_sequences[new_name] = seq

        # Update the key in the divergence data dictionary
        if self._divergence_data and old_name in self._divergence_data:
             data = self._divergence_data.pop(old_name)
             self._divergence_data[new_name] = data

        # Update the map with the new name for this item ID
        self._treeview_id_to_taxon_name[item] = new_name

        # Update the Treeview cell value for the Taxon column (index 0)
        current_values = list(self.divergence_treeview.item(item, 'values'))
        current_values[0] = new_name
        self.divergence_treeview.item(item, values=current_values)


        self.entry_editor.destroy()
        self.entry_editor = None
        self.entry_item = None
        self.entry_column = None

    def remove_selected_taxon_from_results(self):
         """Removes the selected taxon(s) from the results table and internal data."""
         selected_items = self.divergence_treeview.selection()
         if not selected_items:
              messagebox.showwarning("Warning", "Please select one or more taxons in the results table to remove.")
              return

         # Get the current names using the map associated with each item ID
         taxons_to_remove = [self._treeview_id_to_taxon_name.get(item, self.divergence_treeview.item(item, 'values')[0]) for item in selected_items]

         if not taxons_to_remove: return

         confirm = messagebox.askyesno("Confirm Removal", f"Are you sure you want to remove the selected taxon(s) from the results?\n\n" + "\n".join(taxons_to_remove) + "\n\nNote: This does NOT re-run concatenation or update statistics. It only removes the row from the results and internal sequence/divergence data.")
         if not confirm:
              return

         # Process removals - sort selected items in reverse order to avoid issues with deleting
         for item in sorted(selected_items, reverse=True):
              # Get the current name using the map
              taxon_name = self._treeview_id_to_taxon_name.get(item, self.divergence_treeview.item(item, 'values')[0])

              if taxon_name:
                  # Remove from concatenated sequences
                  if self._concatenated_sequences and taxon_name in self._concatenated_sequences:
                       del self._concatenated_sequences[taxon_name]
                  # Remove from divergence data
                  if self._divergence_data and taxon_name in self._divergence_data:
                       del self._divergence_data[taxon_name]
                  # Remove from the Treeview ID map
                  if item in self._treeview_id_to_taxon_name:
                       del self._treeview_id_to_taxon_name[item]

              # Remove from the Treeview display
              self.divergence_treeview.delete(item)

         # Note: Statistics and Partition are NOT updated after taxon removal.
         # A message could be added here to inform the user of this limitation.
         # messagebox.showinfo("Removal Complete", "Selected taxons removed from results display and internal data. Statistics and Partition were not re-calculated.")


    def excise_entire_taxon_from_results(self):
        """Excises (removes) the selected taxon(s) from the results table and internal data."""
        self.remove_selected_taxon_from_results()


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
        self._treeview_id_to_taxon_name = {} # Clear the ID map


        self._concatenated_sequences = None
        self._partition_data = None
        self._statistics = None
        self._divergence_data = None


    def on_reset(self):
        """Resets the input (genes list) and clears results."""
        if self.gene_names or self.raw_gene_contents or self._concatenated_sequences is not None:
            confirm = messagebox.askyesno("Confirm Reset", "Are you sure you want to remove all loaded genes and clear results?")
            if not confirm:
                return

        self.gene_names = []
        self.raw_gene_contents = []
        self.gene_listbox.delete(0, tk.END)

        self.clear_results_display()

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
                    f.write(f">{taxon}\n")
                    wrapped_seq = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
                    f.write(f"{wrapped_seq}\n")

            messagebox.showinfo("Export Successful", f"Concatenated sequences exported to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save FASTA file:\n{e}")
            traceback.print_exc()


    def export_partition(self):
        """Exports the partition data to a NEXUS formatted file (partition block only)."""
        if not self._partition_data or not self._concatenated_sequences:
            messagebox.showwarning("Warning", "No partition data or concatenated sequences available. Run 'Concatenate & Analyze' first.")
            return

        num_taxa = len(self._concatenated_sequences)
        seq_len = len(list(self._concatenated_sequences.values())[0]) if self._concatenated_sequences and list(self._concatenated_sequences.values()) else 0

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
                     if gene_type.upper() == 'DNA':
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          try:
                               start, end = map(int, coord_range.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    f.write(f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n") # Corrected typo: use clean_gene_name
                                    dna_genes_with_partitions_in_partition.append(clean_gene_name)
                               elif length > 0:
                                    f.write(f"  [ Warning: Gene {clean_gene_name} ({coord_range}) is too short ({length} bp) for full codon partitions. ]\n")
                          except ValueError:
                               f.write(f"  [ Warning: Could not parse range '{coord_range}' for codon partitions for gene {clean_gene_name}. ]\n")
                     elif gene_type.upper() == 'PROTEIN':
                          f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")

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

                if dna_genes_with_partitions_in_partition:
                    f.write("\n  [ Link block (by codon position, DNA genes only) ]\n")
                    f.write("  link characters = ")
                    pos1_links = []
                    pos2_links = []
                    pos3_links = []

                    for name, range_str, gene_type in self._partition_data:
                         if gene_type.upper() == 'DNA':
                              clean_name = re.sub(r'\W+', '_', name)
                              try:
                                   start, end = map(int, range_str.split('-'))
                                   length = end - start + 1
                                   if length >= 3:
                                        pos1_links.append(f"{clean_name}_pos1 : {start}-{end}\\3")
                                        pos2_links.append(f"{clean_name}_pos2 : {start+1}-{end}\\3")
                                        pos3_links.append(f"{clean_name}_pos3 : {start+2}-{end}\\3")
                              except ValueError:
                                   pass

                    all_codon_links = pos1_links + pos2_links + pos3_links
                    if all_codon_links:
                        f.write(", ".join(all_codon_links) + ";\n")
                    else:
                         f.write("[ No DNA genes long enough for codon link block ]\n")

                f.write("END; [PAUP]\n\n")

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
                     if gene_type.upper() == 'DNA':
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          try:
                               start, end = map(int, coord_range.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    f.write(f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n") # Corrected typo: use clean_gene_name
                                    f.write(f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n") # Corrected typo: use clean_gene_name
                                    mr_bayes_dna_charsets_exist_mb = True
                                    pos1_names_mb.append(f"{clean_gene_name}_pos1")
                                    pos2_names_mb.append(f"{clean_gene_name}_pos2")
                                    pos3_names_mb.append(f"{clean_gene_name}_pos3")
                          except ValueError:
                               pass
                     elif gene_type.upper() == 'PROTEIN':
                        self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")

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
                          f.write("  partition by_codon_pos = {0}: {1};\n".format(
                               len(partition_groups_mrbayes_mb), ", ".join(partition_groups_mrbayes_mb))
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
                    f.write("  [ No partitions generated ]\n")

                f.write("END; [mrbayes]\n")


            messagebox.showinfo("Export Successful", f"Partition data exported to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save partition file:\n{e}")
            traceback.print_exc()


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
                seq_len = len(list(self._concatenated_sequences.values())[0]) if self._concatenated_sequences and list(self._concatenated_sequences.values()) else 0

                if num_taxa == 0 or seq_len == 0:
                     messagebox.showwarning("Warning", "No concatenated sequences or zero length. Cannot export NEXUS matrix.")
                     return

                f.write("BEGIN TAXA;\n")
                f.write(f"  DIMENSIONS NTAX={num_taxa};\n")
                f.write("  TAXLABELS\n")

                for taxon in sorted_taxons:
                    if re.search(r'[\s\'"`=;:,\[\]\(\)]+', str(taxon)):
                         f.write(f"    '{taxon}'\n")
                    else:
                         f.write(f"    {taxon}\n")
                f.write("  ;\n")
                f.write("END; [TAXA]\n\n")


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
                    taxon_label = f"'{taxon}'" if re.search(r'[\s\'"`=;:,\[\]\(\)]+', str(taxon)) else taxon

                    f.write(f"{taxon_label}\n")
                    wrapped_seq = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
                    indented_wrapped_seq = '\n'.join([f"  {line}" for line in wrapped_seq.splitlines()])
                    f.write(f"{indented_wrapped_seq}\n")


                f.write("  ;\n")
                f.write("END; [DATA]\n\n")

                f.write("BEGIN PAUP;\n")
                f.write("  [ Charsets for gene partitions ]\n")
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                    clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                    f.write(f"  charset {clean_gene_name} = {coord_range};\n")

                f.write("\n  [ Charsets for codon positions (DNA genes only) ]\n")
                dna_genes_with_partitions_in_nexus = []
                for backend_gene_name, coord_range, gene_type in self._partition_data:
                     if gene_type.upper() == 'DNA':
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          try:
                               start, end = map(int, coord_range.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    f.write(f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n") # Corrected typo: use clean_gene_name
                                    dna_genes_with_partitions_in_nexus.append(clean_gene_name)
                               elif length > 0:
                                    f.write(f"  [ Warning: Gene {clean_gene_name} ({coord_range}) is too short ({length} bp) for full codon partitions. ]\n")
                          except ValueError:
                               f.write(f"  [ Warning: Could not parse range '{coord_range}' for codon partitions for gene {clean_gene_name}. ]\n")
                     elif gene_type.upper() == 'PROTEIN':
                           f.write(f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")

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

                if dna_genes_with_partitions_in_nexus:
                    f.write("\n  [ Link block (by codon position, DNA genes only) ]\n")
                    f.write("  link characters = ")
                    pos1_links = []
                    pos2_links = []
                    pos3_links = []

                    for name, range_str, gene_type in self._partition_data:
                         if gene_type.upper() == 'DNA':
                              clean_name = re.sub(r'\W+', '_', name)
                              try:
                                   start, end = map(int, range_str.split('-'))
                                   length = end - start + 1
                                   if length >= 3:
                                        pos1_links.append(f"{clean_name}_pos1 : {start}-{end}\\3")
                                        pos2_links.append(f"{clean_name}_pos2 : {start+1}-{end}\\3")
                                        pos3_links.append(f"{clean_name}_pos3 : {start+2}-{end}\\3")
                              except ValueError:
                                   pass

                    all_codon_links = pos1_links + pos2_links + pos3_links
                    if all_codon_links:
                        f.write(", ".join(all_codon_links) + ";\n")
                    else:
                         f.write("[ No DNA genes long enough for codon link block ]\n")

                f.write("END; [PAUP]\n\n")

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
                     if gene_type.upper() == 'DNA':
                          clean_gene_name = re.sub(r'\W+', '_', backend_gene_name)
                          try:
                               start, end = map(int, coord_range.split('-'))
                               length = end - start + 1
                               if length >= 3:
                                    f.write(f"  charset {clean_gene_name}_pos1 = {start}-{end}\\3;\n")
                                    f.write(f"  charset {clean_gene_name}_pos2 = {start+1}-{end}\\3;\n") # Corrected typo: use clean_gene_name
                                    f.write(f"  charset {clean_gene_name}_pos3 = {start+2}-{end}\\3;\n") # Corrected typo: use clean_gene_name
                                    mr_bayes_dna_charsets_exist_mb = True
                                    pos1_names_mb.append(f"{clean_gene_name}_pos1")
                                    pos2_names_mb.append(f"{clean_name}_pos2") # Fix: use clean_name_pos2
                                    pos3_names_mb.append(f"{clean_name}_pos3") # Fix: use clean_name_pos3
                          except ValueError:
                               pass
                     elif gene_type.upper() == 'PROTEIN':
                        self.partition_text.insert(tk.END, f"  [ Note: Gene {backend_gene_name} is PROTEIN type, no codon partitions generated. ]\n")

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
                          f.write("  partition by_codon_pos = {0}: {1};\n".format(
                               len(partition_groups_mrbayes_mb), ", ".join(partition_groups_mrbayes_mb))
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
                    f.write("  [ No partitions generated ]\n")

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

