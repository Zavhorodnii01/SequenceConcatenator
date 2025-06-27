import re
import math
import sys # Import sys for stderr output
import traceback


class SequenceConcatenator:
    """
    A class to parse gene sequence files (FASTA, Nexus, GenBank),
    concatenate sequences for common taxons, and generate statistics
    and partition information, including per-taxon per-gene divergence.
    """

    # Accept gene_names from the frontend
    def __init__(self, gene_file_contents: list[list[str]], gene_names: list[str]):
        """
        Initializes the SequenceConcatenator with raw gene file contents and their names.

        Args:
            gene_file_contents: A list where each element is a list of strings,
                                representing the lines of a single gene file.
            gene_names: A list of strings, representing the names assigned to each gene file.
                        Must be the same length as gene_file_contents.
        """
        if len(gene_file_contents) != len(gene_names):
             # This check should ideally prevent issues, caught in the frontend
             print("Error: Length of gene_file_contents and gene_names mismatch.", file=sys.stderr)
             # Raise an error as this indicates a serious sync issue
             raise ValueError("Length of gene_file_contents and gene_names must be the same.")

        self.__raw_gene_contents = gene_file_contents
        self.__frontend_gene_names = gene_names # Store the names passed from frontend
        self.__parsed_gene_data = [] # [{'taxon1': 'seq1', 'taxon2': 'seq2', ...}, ...] list of dicts per gene
        self.__gene_info = [] # Stores info like {'name': 'GeneName', 'type': 'DNA', 'length': 100, 'start': 1, 'end': 100}
        self.__all_taxons = [] # All taxons found across all input files (before concatenation filtering)

        self.__concatenated_sequences = {} # {taxon_name: concatenated_sequence, ...} - Only taxons with data after concat
        self.__partition_data = [] # [(gene_name, 'start-end', gene_type), ...]
        self.__statistics = {} # Dictionary holding various stats, including initial divergence


        self._parse_all_genes()
        self._collect_all_taxons() # Collects all taxons from parsed data
        # _concatenate_sequences filters to taxons with data and sets lengths/positions in __gene_info
        self.__concatenated_sequences = self._concatenate_sequences()

        # Filter __all_taxons to only include those that actually ended up in concatenated sequences
        self.__all_taxons = sorted(list(self.__concatenated_sequences.keys()))

        # Only keep gene_info entries that resulted in a non-zero length in the concatenation
        self.__gene_info = [info for info in self.__gene_info if info.get('length', 0) > 0]


        # _calculate_partition uses the updated __gene_info
        self.__partition_data = self._calculate_partition()
        # _calculate_statistics calculates various stats including the *initial* divergence relative to the first taxon
        self.__statistics = self._calculate_statistics() # This call now uses the modified _perform_divergence_calculation which defaults to the first taxon


    def get_concatenated_sequences(self) -> dict[str, str]:
        """
        Returns the dictionary of concatenated sequences.
        """
        return self.__concatenated_sequences

    def get_statistics(self) -> dict:
        """
        Returns the calculated statistics for the concatenated alignment.
        """
        return self.__statistics

    def get_partition(self) -> list[tuple[str, str, str]]:
        """
        Returns the partition data for the concatenated alignment.
        """
        return self.__partition_data

    # Added a public method to allow the UI to get the list of genes with length > 0
    def get_processed_gene_info(self) -> list[dict]:
         """
         Returns the list of gene info dictionaries for genes included in the concatenated alignment.
         """
         # Return a copy to prevent external modification
         return list(self.__gene_info)


    # New public method to recalculate divergence using the instance's *internal* data
    def recalculate_divergence_using_internal_data(self, reference_taxon_name: str = None) -> dict:
        """
        Calculates divergence using the instance's concatenated sequences and gene info,
        relative to the specified reference taxon.

        Args:
            reference_taxon_name: The name of the taxon to use as the reference.
                                If None or not found in sequences, the first taxon is used.

        Returns:
            A dictionary containing recalculatd divergence data per taxon.
        """
        # Use the internal data and pass them to the core calculation logic
        # Use the final list of taxons from concatenated sequences
        taxons_in_concat = sorted(list(self.__concatenated_sequences.keys()))

        return self._perform_divergence_calculation(
            self.__concatenated_sequences,
            self.__gene_info, # Pass the processed gene_info (only genes with length > 0)
            taxons_in_concat, # Pass current taxons in concat
            reference_taxon_name
        )


    # Internal parsing methods remain the same
    def _parse_all_genes(self) -> None:
        """
        Parses all raw gene file contents, detects format, and stores data and basic info.
        Uses the gene name provided by the frontend.
        Initial __gene_info list is built here, to be refined by _concatenate_sequences.
        """
        self.__parsed_gene_data = [] # Clear existing data on re-parse
        self.__gene_info = []      # Clear existing info on re-parse

        for i, file_lines in enumerate(self.__raw_gene_contents):
            content_str = "".join(file_lines).strip()
            current_gene_name = self.__frontend_gene_names[i] if i < len(self.__frontend_gene_names) else f"gene_index_{i+1}"

            parsed_data = {}
            determined_seq_type = 'Unknown'

            try:
                is_nexus = content_str.lower().strip().startswith('#nexus') or 'begin data' in content_str.lower()
                is_genbank = content_str.lower().strip().startswith('locus') or 'origin' in content_str.lower() or 'version' in content_str.lower() or 'ncbi' in content_str.lower()
                is_fasta = content_str.strip().startswith('>')

                if is_nexus:
                    parsed_data = self._parse_nexus(file_lines)
                elif is_genbank:
                    parsed_data = self._parse_genbank(file_lines)
                elif is_fasta:
                    parsed_data = self._parse_fasta(file_lines)
                else:
                     # Try FASTA as a default
                    parsed_data = self._parse_fasta(file_lines)
                    if not parsed_data:
                         # Simple single-sequence fallback
                         non_empty_lines = [line.strip() for line in file_lines if line.strip() and not line.strip().startswith('#') and not line.strip().startswith('[')]
                         if non_empty_lines:
                              seq_content = "".join(non_empty_lines).replace(' ', '').replace('\t','')
                              if seq_content:
                                  dummy_taxon_name = f"Taxon_for_{current_gene_name}_File"
                                  print(f"Info: Gene file index {i+1} ('{current_gene_name}') contains data but no recognized header. Treating as simple single-sequence for taxon '{dummy_taxon_name}'.", file=sys.stderr)
                                  parsed_data = {dummy_taxon_name: seq_content}
                                  if not seq_content:
                                      print(f"Warning: Gene file index {i+1} ('{current_gene_name}') had lines but yielded no sequence data after cleaning.", file=sys.stderr)
                                      parsed_data = {}

                # Determine sequence type
                if parsed_data:
                    first_seq = next((re.sub(r'[\s\-?]+', '', seq).upper() for seq in parsed_data.values() if seq), "")

                    if first_seq:
                         dna_chars = set("ACGTU")
                         protein_chars = set("ACDEFGHIKLMNPQRSTVWY")
                         ambiguous_dna = set("RYSWKMBDHVN")
                         ambiguous_protein = set("BJZX")
                         all_dna_chars_strict = dna_chars | ambiguous_dna
                         all_protein_chars_strict = protein_chars | ambiguous_protein
                         seq_chars_present = set(first_seq)

                         is_potential_dna = seq_chars_present.issubset(all_dna_chars_strict)
                         is_potential_protein = seq_chars_present.issubset(all_protein_chars_strict)

                         if not seq_chars_present: determined_seq_type = 'Unknown'
                         elif is_potential_dna and not (seq_chars_present & protein_chars): determined_seq_type = 'DNA'
                         elif is_potential_protein and not (seq_chars_present & dna_chars): determined_seq_type = 'Protein'
                         elif is_potential_dna and is_potential_protein:
                              if seq_chars_present & set("TU"): determined_seq_type = 'DNA'
                              elif seq_chars_present & set("FILPQEKRWYV"): determined_seq_type = 'Protein'
                              else: determined_seq_type = 'Unknown'
                         else: determined_seq_type = 'Unknown'

                    if parsed_data:
                         self.__parsed_gene_data.append(parsed_data)
                         # Initialize gene info with length 0 and placeholder positions
                         self.__gene_info.append({'name': current_gene_name, 'type': determined_seq_type, 'length': 0, 'start': 0, 'end': 0})
                    else:
                         print(f"Warning: Gene file index {i+1} ('{current_gene_name}') parsed but yielded no sequences.", file=sys.stderr)

            except Exception as e:
                print(f"Error parsing gene file index {i+1} ('{current_gene_name}'): {e}", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)

        # Filter out genes that resulted in no parsed data dictionaries
        valid_genes_data = []
        initial_gene_info_before_concat = []
        # We need to map the index of the parsed data dict back to the original gene_info
        # to get the name and type.
        original_gene_info_map = {i: {'name': self.__gene_info[i].get('name', f'UnnamedGene_{i+1}'),
                                      'type': self.__gene_info[i].get('type', 'Unknown')}
                                  for i in range(len(self.__gene_info))}

        for i, gene_data in enumerate(self.__parsed_gene_data):
             if gene_data:
                  valid_genes_data.append(gene_data)
                  # Use the original info based on index
                  initial_gene_info_before_concat.append(original_gene_info_map.get(i, {'name': f'UnnamedGene_{i+1}_Error', 'type': 'Unknown'}))


        self.__parsed_gene_data = valid_genes_data
        self.__gene_info = initial_gene_info_before_concat # Use this list to be updated in concat


    def _parse_fasta(self, lines: list[str]) -> dict[str, str]:
        sequences = {}
        current_name = None
        current_seq_lines = []
        processed_lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]

        for line in processed_lines:
            if line.startswith(">"):
                if current_name is not None and current_seq_lines:
                    sequences[current_name] = "".join(current_seq_lines)
                match = re.match(r'>\s*(.+)', line)
                if match: current_name = match.group(1).strip()
                else: current_name = line[1:].strip()
                if not current_name:
                     current_name = f"UnnamedTaxon_{len(sequences) + 1}"
                     print(f"Warning: Found FASTA header with no name. Assigning default name '{current_name}'.", file=sys.stderr)
                current_seq_lines = []
            elif current_name is not None:
                clean_seq_line = re.sub(r'[\s\d]+', '', line)
                current_seq_lines.append(clean_seq_line)
        if current_name is not None and current_seq_lines:
            sequences[current_name] = "".join(current_seq_lines)
        return sequences

    def _parse_nexus(self, lines: list[str]) -> dict[str, str]:
        sequences = {}
        in_matrix = False
        taxon_sequences_buffer = {}
        current_taxon = None
        inside_block_level = 0

        for raw_line in lines:
            line = raw_line.strip()
            if not line: continue
            if '[' in line and ']' in line:
                 line = re.sub(r'\[.*?\]', '', line).strip()
                 if not line: continue
            line_lower = line.lower()

            if line_lower.startswith('begin'):
                 inside_block_level += 1
                 if ' data;' in line_lower and inside_block_level == 1:
                      sequences = {}
                      in_matrix = False
                      taxon_sequences_buffer = {}
                      current_taxon = None
                 continue
            if line_lower.startswith('end;'):
                 if in_matrix:
                      if current_taxon is not None and current_taxon in taxon_sequences_buffer:
                           sequences[current_taxon] = "".join(taxon_sequences_buffer[current_taxon])
                           del taxon_sequences_buffer[current_taxon]
                      for t, seq_parts in list(taxon_sequences_buffer.items()):
                           if t not in sequences:
                              sequences[t] = "".join(seq_parts)
                           del taxon_sequences_buffer[t]
                      in_matrix = False
                 inside_block_level = max(0, inside_block_level - 1)
                 continue
            if in_matrix:
                match = re.match(r"^\s*(?:(['\"])(.*?)\1|([^'\s]+))\s*(\S+)\s*;", line)
                if not match: match = re.match(r"^\s*(?:(['\"])(.*?)\1|([^'\s]+))\s*(\S+)", line)
                if match:
                    quoted_name = match.group(2)
                    unquoted_name = match.group(3)
                    taxon_name = quoted_name if quoted_name is not None else unquoted_name
                    sequence_part = match.group(4).strip()
                    if taxon_name:
                        if current_taxon is not None and current_taxon != taxon_name and current_taxon in taxon_sequences_buffer:
                             if taxon_sequences_buffer[current_taxon]:
                                  if current_taxon not in sequences:
                                       sequences[current_taxon] = "".join(taxon_sequences_buffer[current_taxon])
                        current_taxon = taxon_name
                        if taxon_name not in taxon_sequences_buffer: taxon_sequences_buffer[taxon_name] = []
                        taxon_sequences_buffer[taxon_name].append(sequence_part)
                elif current_taxon is not None:
                     sequence_part = line.replace(';', '').strip()
                     if sequence_part:
                          if current_taxon in taxon_sequences_buffer:
                               taxon_sequences_buffer[current_taxon].append(sequence_part)

        for t, seq_parts in taxon_sequences_buffer.items():
             if t not in sequences:
                  sequences[t] = "".join(seq_parts)

        cleaned_sequences = {name: re.sub(r'[\s\d]+', '', seq).upper() for name, seq in sequences.items() if seq.strip()}
        return cleaned_sequences


    def _parse_genbank(self, lines: list[str]) -> dict[str, str]:
        sequences = {}
        current_organism_name = None
        in_origin = False
        sequence_parts = []

        for line in lines:
            line_stripped = line.rstrip()
            if not line_stripped.strip(): continue

            if line_stripped.strip() == '//':
                if current_organism_name is not None and sequence_parts:
                    sequence = "".join(sequence_parts)
                    clean_sequence = re.sub(r'[\s\d]+', '', sequence).upper()
                    if clean_sequence:
                         sequences[current_organism_name] = clean_sequence
                current_organism_name = None
                in_origin = False
                sequence_parts = []
                continue

            org_match = re.match(r'^\s{2,}ORGANISM\s+(.*?)\s*\.?$', line_stripped, re.IGNORECASE)
            if org_match:
                current_organism_name = org_match.group(1).strip()
                in_origin = False
                sequence_parts = []
                continue

            if line_stripped.strip().lower() == 'origin':
                in_origin = True
                sequence_parts = []
                continue

            if in_origin:
                sequence_parts.append(line_stripped)
                continue

            if current_organism_name is None:
                 org_feature_match = re.search(r'/organism="([^"]+)"', line_stripped, re.IGNORECASE)
                 if org_feature_match:
                      current_organism_name = org_feature_match.group(1).strip()

        if current_organism_name is not None and sequence_parts:
             sequence = "".join(sequence_parts)
             clean_sequence = re.sub(r'[\s\d]+', '', sequence).upper()
             if clean_sequence:
                   sequences[current_organism_name] = clean_sequence

        cleaned_sequences = {name: seq for name, seq in sequences.items() if seq}
        return cleaned_sequences


    def _collect_all_taxons(self) -> None:
        """
        Collects all unique taxon names from the parsed data dictionaries.
        This list includes all taxons found in *any* gene file *that was successfully parsed*.
        This list is used as the basis for creating the concatenated sequences dictionary keys.
        """
        all_taxons_set = set()
        # Iterate through the data dictionaries that were successfully parsed
        for gene_data_dict in self.__parsed_gene_data:
            all_taxons_set.update(gene_data_dict.keys())
        self.__all_taxons = sorted(list(all_taxons_set))


    def _concatenate_sequences(self) -> dict[str, str]:
        """
        Concatenates sequences for each taxon found in __all_taxons.
        Adds gaps for missing data or sequences of unexpected length.
        Updates the length and position info in __gene_info for genes
        that successfully contributed sequence data.
        Filters the final concatenated sequences to only include taxons that were in __all_taxons
         and had non-empty sequences after concatenation.
        """
        # Initialize concatenated sequences for all taxons collected from *successfully parsed* genes
        concatenated_sequences = {taxon: "" for taxon in self.__all_taxons}
        current_position = 1 # 1-based indexing for partition/gene_info

        # We will build a *new* gene_info list containing only genes that had length > 0
        # and updating their start/end positions in the concatenated alignment.
        updated_gene_info = []
        # Create a mapping from original gene_info index to gene name/type
        # Use a copy of the original __gene_info list
        original_gene_info_list_copy = list(self.__gene_info)


        # Iterate through the parsed data dictionaries (only those that were successfully parsed)
        for i, gene_data_dict in enumerate(self.__parsed_gene_data):
             # Get the name and type from the original gene_info list using the index
             # This assumes __parsed_gene_data and the original_gene_info_list_copy are aligned by index
             # We should verify the index is valid
             if i >= len(original_gene_info_list_copy):
                  print(f"Internal Error: Mismatch indexing original gene info list at index {i}.", file=sys.stderr)
                  # Fallback name/type
                  gene_name = f'UnnamedGene_{i+1}_Error'
                  gene_type = 'Unknown'
             else:
                  gene_info_entry = original_gene_info_list_copy[i] # Use the original info
                  gene_name = gene_info_entry.get('name', f'UnnamedGene_{i+1}')
                  gene_type = gene_info_entry.get('type', 'Unknown')


             # Determine the *actual* length of this gene segment from the parsed data
             gene_length = 0
             if gene_data_dict:
                  first_non_empty_seq = next((seq for seq in gene_data_dict.values() if seq), "")
                  gene_length = len(first_non_empty_seq)
                  # Note: Assumes all sequences for a given gene in the source file should have the same length.


             if gene_length > 0:
                  # This gene contributes to the concatenated alignment.
                  # Update the gene info with the determined length and position range
                  gene_start = current_position
                  gene_end = current_position + gene_length - 1
                  updated_gene_info.append({'name': gene_name, 'type': gene_type, 'length': gene_length, 'start': gene_start, 'end': gene_end})

                  # Now, append the sequence (or gaps) for each taxon to the concatenated sequence
                  for taxon in self.__all_taxons: # Iterate through ALL taxons found initially from *parsed* data
                      sequence = gene_data_dict.get(taxon) # Get sequence for this taxon from this gene's data

                      if sequence is not None:
                          # If taxon was present in this gene's data dictionary
                          if len(sequence) == gene_length:
                              concatenated_sequences[taxon] += sequence
                          else:
                              # Length mismatch - pad or truncate
                              print(f"Warning: Taxon '{taxon}' sequence for gene '{gene_name}' has unexpected length ({len(sequence)} vs {gene_length}). Padding/Truncating.", file=sys.stderr)
                              if len(sequence) < gene_length:
                                  concatenated_sequences[taxon] += sequence + "-" * (gene_length - len(sequence))
                              else: # len(sequence) > gene_length
                                  concatenated_sequences[taxon] += sequence[:gene_length]
                      else:
                          # Taxon is missing from this gene_data_dict, append gaps of the expected length
                          concatenated_sequences[taxon] += "-" * gene_length

                  # Update the current position for the next gene ONLY if this gene had non-zero length
                  current_position += gene_length

             else:
                  print(f"Info: Gene '{gene_name}' has effective length 0 after parsing or contains no sequences. Not added to concatenated alignment or partition.", file=sys.stderr)


        # Update the official __gene_info list to contain only genes with length > 0 and updated positions
        self.__gene_info = updated_gene_info

        # Filter out taxons that ended up with zero length sequences (e.g., if all genes failed to parse for that taxon)
        # Or simply keep taxons whose sequence length matches the total alignment length (due to padding)
        expected_total_length = sum(info.get('length', 0) for info in self.__gene_info)
        final_concatenated_sequences = {
            taxon: seq for taxon, seq in concatenated_sequences.items() if len(seq) == expected_total_length # Only keep if sequence length matches
        }

        # Update __all_taxons to reflect only taxons actually present in the final alignment
        self.__all_taxons = sorted(list(final_concatenated_sequences.keys()))

        # Final check: ensure all remaining concatenated sequences have the same total length (should be true due to filtering)
        if self.__all_taxons:
             sample_seq_len = len(next(iter(final_concatenated_sequences.values())))
             if sample_seq_len != expected_total_length:
                  print(f"Internal Error: Final concatenated sequences have inconsistent total length ({sample_seq_len} vs expected {expected_total_length}). This indicates a bug in concatenation logic.", file=sys.stderr)


        return final_concatenated_sequences


    def _calculate_partition(self) -> list[tuple[str, str, str]]:
        """
        Generates partition information based on the final __gene_info list
        (containing only genes that contributed length > 0).
        Positions are 1-based.
        """
        partition = []
        # Use the __gene_info list updated by _concatenate_sequences
        # This list already contains only genes with length > 0 and valid ranges
        for gene_info in self.__gene_info:
            gene_name = gene_info['name']
            start = gene_info['start']
            end = gene_info['end']
            gene_type = gene_info.get('type', 'Unknown')

            partition.append((gene_name, f"{start}-{end}", gene_type))

        return partition


    def _calculate_statistics(self) -> dict:
        """
        Calculates various statistics, including per-taxon per-gene divergence.
        This method uses the final concatenated sequences and gene info.
        Initial divergence is calculated relative to the first taxon in __all_taxons.
        """
        statistics = {}

        # Calculate statistics based on the FINAL concatenated sequences and gene info
        taxons_in_concat = self.__all_taxons # Use the final list of taxons from concatenated sequences
        num_taxa = len(taxons_in_concat)
        statistics["Number of Taxa"] = num_taxa
        # Count genes that contributed sequence (those in the final __gene_info list)
        statistics["Number of Genes"] = len(self.__gene_info)


        if num_taxa == 0 or not self.__concatenated_sequences:
             statistics["Total Length"] = 0
             statistics["Percentage Overall Missing Data (%)"] = 0.0
             statistics["Missing Data per Taxon (count)"] = {}
             statistics["Missing Data per Gene (count)"] = {}
             statistics["Gene Lengths (in concatenated alignment)"] = {}
             statistics["Taxon-Gene Matrix Sparsity (%)"] = 0.0
             statistics["Divergence Data"] = {} # Initial divergence data is empty if no taxons/data
             return statistics

        total_length = len(next(iter(self.__concatenated_sequences.values()))) # Safely get length
        statistics["Total Length"] = total_length

        total_cells = num_taxa * total_length if total_length > 0 else 0
        total_missing_chars = 0
        missing_data_per_taxon = {}

        for taxon in taxons_in_concat:
            seq = self.__concatenated_sequences.get(taxon, "")
            missing_count = seq.count("-") + seq.count("?")
            missing_data_per_taxon[taxon] = missing_count
            total_missing_chars += missing_count

        statistics["Missing Data per Taxon (count)"] = missing_data_per_taxon


        # Calculate missing data per gene based on concatenated sequences
        # Use the final __gene_info list (only genes with length > 0)
        missing_data_per_gene = {info['name']: 0 for info in self.__gene_info}
        for gene_info in self.__gene_info:
             gene_name = gene_info['name']
             start_0based = gene_info['start'] - 1
             end_0based = gene_info['end']
             # gene_segment_length = gene_info['length'] # No need, use slice length


             if start_0based >= 0 and end_0based >= start_0based and total_length >= end_0based:
                  for taxon in taxons_in_concat:
                       seq = self.__concatenated_sequences.get(taxon, "")
                       if seq and len(seq) >= end_0based:
                            seq_segment = seq[start_0based:end_0based]
                            missing_count_segment = seq_segment.count("-") + seq_segment.count("?")
                            missing_data_per_gene[gene_name] += missing_count_segment
                       # else: Warning already printed in concatenation if padding failed

        statistics["Missing Data per Gene (count)"] = missing_data_per_gene

        if total_cells > 0:
            percentage_missing_overall = (total_missing_chars / total_cells) * 100
            statistics["Percentage Overall Missing Data (%)"] = round(percentage_missing_overall, 2)
        else:
            statistics["Percentage Overall Missing Data (%)"] = 0.0

        # Calculate Taxon-Gene Matrix Sparsity
        num_genes_in_concat = len(self.__gene_info) # Number of genes with length > 0
        num_expected_gene_segments = num_taxa * num_genes_in_concat if num_genes_in_concat > 0 else 0
        num_fully_missing_gene_segments = 0

        if num_expected_gene_segments > 0:
             for taxon in taxons_in_concat:
                  seq = self.__concatenated_sequences.get(taxon, "")
                  for gene_info in self.__gene_info:
                       gene_segment_length = gene_info['length']
                       start_0based = gene_info['start'] - 1
                       end_0based = gene_info['end']

                       if seq and len(seq) >= end_0based:
                            gene_segment = seq[start_0based:end_0based]
                            if len(gene_segment) == gene_segment_length and all(char in ('-', '?') for char in gene_segment):
                                 num_fully_missing_gene_segments += 1
                       # else: Warning already printed in concatenation if padding failed

        if num_expected_gene_segments > 0:
             sparsity_percentage = (num_fully_missing_gene_segments / num_expected_gene_segments) * 100
             statistics["Taxon-Gene Matrix Sparsity (%)"] = round(sparsity_percentage, 2)
        else:
             statistics["Taxon-Gene Matrix Sparsity (%)"] = 0.0


        # Gene lengths are already in __gene_info, just format for stats output
        gene_lengths_in_concat = {info['name']: info['length'] for info in self.__gene_info}
        statistics["Gene Lengths (in concatenated alignment)"] = gene_lengths_in_concat


        # Calculate the *initial* divergence data relative to the first taxon
        initial_reference_taxon = taxons_in_concat[0] if taxons_in_concat else None
        if initial_reference_taxon:
            # Use the core calculation logic
            divergence_data = self._perform_divergence_calculation(
                self.__concatenated_sequences,
                self.__gene_info, # Pass the processed gene_info
                taxons_in_concat, # Pass current taxons in concat
                initial_reference_taxon # Calculate relative to the first taxon
            )
            statistics["Divergence Data"] = divergence_data
        else:
             statistics["Divergence Data"] = {}


        return statistics

    # Renamed core calculation logic and made it accept all necessary data as arguments
    # This method is now called by _calculate_statistics (for initial calculation)
    # and by recalculate_divergence_using_internal_data (for subsequent calculations)
    def _perform_divergence_calculation(self, concatenated_sequences: dict[str, str], gene_info: list[dict], taxons_in_concat: list[str], reference_taxon_name: str = None) -> dict:
        """
        Core logic to calculate difference statistics for each taxon relative to a specified
        reference taxon for each gene segment and overall.

        Args:
            concatenated_sequences: The dictionary of concatenated sequences.
            gene_info: The list of gene info dictionaries (name, start, end, type, length) for genes > 0 length.
            taxons_in_concat: The list of taxons present in the concatenated sequences.
            reference_taxon_name: The name of the taxon to use as the reference.
                                If None or not found in sequences, the first taxon in taxons_in_concat is used.

        Returns:
            A dictionary containing divergence data per taxon.
        """
        divergence_data = {}

        # Use the list of taxons actually present in the concatenated alignment
        taxons = taxons_in_concat

        if not taxons or len(taxons) < 1 or not concatenated_sequences or not gene_info:
            # If no taxons or no data, return empty dict or dicts initialized for taxons
            for taxon in taxons: # Will iterate only if taxons is not empty
                 divergence_data[taxon] = {
                     'Total score': 0.0,
                     'No of charsets': 0,
                 }
            return divergence_data

        # Determine the actual reference taxon to use
        if reference_taxon_name is not None and reference_taxon_name in taxons:
             reference_taxon = reference_taxon_name
        else:
             reference_taxon = taxons[0] # Default to the first taxon if specified is invalid or None

        ref_seq = concatenated_sequences.get(reference_taxon, "")
        # total_alignment_length = len(ref_seq) # Not needed directly in this function


        # Initialize data structure for each taxon
        for taxon in taxons:
             divergence_data[taxon] = {
                  'Total score': 0.0,
                  'No of charsets': 0, # Counts genes where this taxon has non-gap data
             }
             # Initialize gene scores as "N/A" using the actual gene names from gene_info
             for gene in gene_info: # gene_info already only includes genes with length > 0
                 divergence_data[taxon][gene['name']] = "N/A"


        # Calculate per-gene differences and count non-gap charsets
        genes_for_divergence = gene_info # Use the filtered list


        for gene in genes_for_divergence:
             gene_name = gene['name']
             start_0based = gene['start'] - 1
             end_0based = gene['end']
             gene_segment_length = gene['length'] # Use length from gene_info

             # Get reference segment
             ref_seq_full = concatenated_sequences.get(reference_taxon, "")
             # This check should ideally pass if data comes from a valid concatenated_sequences dict
             if len(ref_seq_full) < end_0based:
                  print(f"Internal Error: Reference taxon '{reference_taxon}' sequence unexpectedly shorter ({len(ref_seq_full)}) than end position ({end_0based}) for gene '{gene_name}'. Cannot calculate divergence for this gene.", file=sys.stderr)
                  for taxon in taxons: # Mark this gene as error for all taxons
                       divergence_data[taxon][gene_name] = "Error (Ref Seq Short)"
                  continue
             ref_segment = ref_seq_full[start_0based:end_0based]


             # Calculate scores for each taxon relative to the reference taxon for this gene segment
             for taxon in taxons:
                  # Reference taxon score is always 0 for its gene column
                  if taxon == reference_taxon:
                       # For gene columns in the reference row, it's 0%
                       divergence_data[taxon][gene_name] = f"0% #0 (0.00%)"
                       continue # Skip detailed comparison for reference taxon vs itself


                  taxon_seq = concatenated_sequences.get(taxon, "")
                  if len(taxon_seq) < end_0based:
                       print(f"Internal Error: Taxon '{taxon}' sequence unexpectedly shorter ({len(taxon_seq)}) than end position ({end_0based}) for gene '{gene_name}'. Cannot calculate divergence for this gene.", file=sys.stderr)
                       divergence_data[taxon][gene_name] = "Error (Seq Short)"
                       continue
                  taxon_segment = taxon_seq[start_0based:end_0based]

                  diff_count = 0
                  taxon_has_non_gap_in_gene_segment = False

                  for j in range(gene_segment_length):
                       ref_char = ref_segment[j]
                       taxon_char = taxon_segment[j]

                       if taxon_char not in ('-', '?'):
                            taxon_has_non_gap_in_gene_segment = True

                       if ref_char not in ('-', '?') and taxon_char not in ('-', '?'):
                            if ref_char.upper() != taxon_char.upper():
                                 diff_count += 1

                  # Calculate percentage based on gene segment length
                  percentage = (diff_count / gene_segment_length) * 100 if gene_segment_length > 0 else 0.0

                  leading_percent_int = int(round(percentage))
                  leading_percent_int = max(0, min(99, leading_percent_int)) # Cap 0-99

                  divergence_data[taxon][gene_name] = f"{leading_percent_int}% #{diff_count} ({percentage:.2f}%)"

                  # 'No of charsets' count is done below


        # Calculate Total Score per taxon (Average of gene percentages where taxon had data)
        # And finalize the 'No of charsets' count per taxon
        genes_that_contributed_length = gene_info # This list already has genes with length > 0
        num_genes_in_calc = len(genes_that_contributed_length)

        for taxon in taxons:
             total_percentage_sum = 0.0
             genes_counted_for_total_score = 0
             taxon_genes_with_data_count = 0 # Recalculate No of charsets here

             for gene in genes_that_contributed_length: # Iterate through genes that actually had length > 0
                  gene_name = gene['name']
                  gene_score_str = divergence_data[taxon].get(gene_name)

                  # Check if the score was successfully calculated and get the precise percentage
                  if gene_score_str and "N/A" not in gene_score_str and "Error" not in gene_score_str:
                       match = re.search(r'\((\d+\.?\d*)%\)', gene_score_str)
                       if match:
                            try:
                                 percentage_value = float(match.group(1))
                                 total_percentage_sum += percentage_value
                                 genes_counted_for_total_score += 1
                            except ValueError:
                                 print(f"Warning: Could not parse percentage from gene score string for taxon '{taxon}', gene '{gene_name}': {gene_score_str}. Skipping for total score calculation.", file=sys.stderr)

                       # Check if this gene contributed non-gap data for the taxon to count charsets
                       taxon_seq = concatenated_sequences.get(taxon, "")
                       start_0based = gene['start'] - 1
                       end_0based = gene['end']
                       # Ensure sequence is long enough before slicing
                       if taxon_seq and len(taxon_seq) >= end_0based:
                            taxon_segment = taxon_seq[start_0based:end_0based]
                            if any(char not in ('-', '?') for char in taxon_segment):
                                 taxon_genes_with_data_count += 1
                       else:
                            # This might happen if a taxon exists in __all_taxons but somehow didn't get a full padded sequence in concat
                            print(f"Internal Warning: Taxon '{taxon}' sequence length unexpected for gene '{gene_name}' during charset count.")


             average_percentage = (total_percentage_sum / genes_counted_for_total_score) if genes_counted_for_total_score > 0 else 0.0
             divergence_data[taxon]['Total score'] = round(average_percentage, 2)

             # Update the 'No of charsets' count for this taxon
             # For the reference taxon, No of charsets is total number of genes with length > 0
             # For others, it's the number of genes where *this taxon* had non-gap data
             if taxon == reference_taxon:
                 divergence_data[taxon]['No of charsets'] = num_genes_in_calc
             else:
                 divergence_data[taxon]['No of charsets'] = taxon_genes_with_data_count


        return divergence_data


# --- Main execution ---
# This section is only for standalone testing of the backend class if needed
if __name__ == "__main__":
    # Example usage - This will only run if you execute SequenceConcatenator.py directly
    print("Running SequenceConcatenator.py as a standalone script (for testing backend logic).")
    print("To run the GUI application, execute SequenceConcatenatorApp.py.")

    # Create dummy raw content (list of lists of strings)
    dummy_fasta_gene1 = [
        ">TaxonA\nACGTACGTACGT", # Length 12
        ">TaxonB\nACGTACGTACGT",
        ">TaxonC\nACGTACGTACGT"
    ] # Gene Name: Gene_One

    dummy_fasta_gene2 = [
        ">TaxonA\nTTTT----GGGG", # Length 12
        ">TaxonB\nCCCCAAAA----",
        ">TaxonD\nGGGGGGGGGGGG" # TaxonD only in gene2
    ] # Gene Name: ImportantGene

    dummy_fasta_gene3 = [
        ">TaxonA\nAAAAAAAAAA", # Length 10
        ">TaxonB\nTTTTTTTTTT",
        ">TaxonC\nCCCCCCCCCC",
        ">TaxonE\nNNNNNNNNNN" # TaxonE only in gene3
    ] # Gene Name: MyGene_v2

     # Gene with no sequences
    dummy_empty_gene = [
        "# This file is empty"
    ] # Gene Name: EmptyGene

    # Raw contents and their assigned names (must match index)
    raw_contents = [
        dummy_fasta_gene1,
        dummy_fasta_gene2,
        dummy_fasta_gene3,
        dummy_empty_gene, # Include an empty one to test filtering
    ]

    assigned_gene_names = [
        "Gene_One",
        "ImportantGene",
        "MyGene_v2",
        "EmptyGene", # Name for the empty gene
    ]

    try:
        # Create backend instance
        concatenator = SequenceConcatenator(raw_contents, assigned_gene_names)

        # Get initial results
        concatenated = concatenator.get_concatenated_sequences()
        stats = concatenator.get_statistics()
        partition = concatenator.get_partition()
        backend_gene_info = concatenator.get_processed_gene_info() # Get the filtered/processed gene info
        initial_divergence = stats.get('Divergence Data', {})

        print("\nConcatenated Sequences:")
        if concatenated:
             for taxon, seq in concatenated.items():
                 print(f">{taxon}\n{seq}")
        else:
             print("No concatenated sequences produced.")

        print("\nStatistics:")
        import json
        class SetEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, set):
                    return list(obj)
                try:
                    return json.JSONEncoder.default(self, obj)
                except TypeError:
                    return str(obj)

        print(json.dumps(stats, indent=2, cls=SetEncoder))

        print("\nPartition:")
        if partition:
             for gene_name, pos_range, gene_type in partition:
                 print(f"Gene: {gene_name}, Range: {pos_range}, Type: {gene_type}")
        else:
             print("No partition data produced.")

        print("\nBackend Processed Gene Info (genes with length > 0):")
        for info in backend_gene_info:
             print(info)

        # Verify divergence data uses correct gene names and try recalculating
        taxons_in_concat = sorted(list(concatenated.keys()))
        if initial_divergence:
             print("\nInitial Divergence Data (Reference: First Taxon):")
             initial_ref = taxons_in_concat[0] if taxons_in_concat else None
             print(f"(Calculated relative to: {initial_ref})")
             print(json.dumps(initial_divergence, indent=2, cls=SetEncoder))


             print("\nRecalculating Divergence Data (Reference: TaxonC):")
             if 'TaxonC' in taxons_in_concat:
                  # Use the public method to recalculate
                  recalculated_divergence = concatenator.recalculate_divergence_using_internal_data('TaxonC')
                  print(json.dumps(recalculated_divergence, indent=2, cls=SetEncoder))

                  print("\nRecalculating Divergence Data (Reference: TaxonD):")
                  if 'TaxonD' in taxons_in_concat:
                       # Use the public method to recalculate
                       recalculated_divergence_d = concatenator.recalculate_divergence_using_internal_data('TaxonD')
                       print(json.dumps(recalculated_divergence_d, indent=2, cls=SetEncoder))
                  else:
                       print("TaxonD not found in concatenated sequences, cannot recalculate relative to TaxonD.")

             else:
                  print("TaxonC not found in concatenated sequences, cannot recalculate relative to TaxonC.")


        else:
             print("\nNo Divergence Data produced.")


    except ValueError as e:
        print(f"Configuration Error: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during backend processing:\n{e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
