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

    def __init__(self, gene_file_contents: list[list[str]]):
        """
        Initializes the SequenceConcatenator with raw gene file contents.

        Args:
            gene_file_contents: A list where each element is a list of strings,
                                representing the lines of a single gene file.
        """
        self.__raw_gene_contents = gene_file_contents
        self.__parsed_gene_data = [] # [{'taxon1': 'seq1', 'taxon2': 'seq2', ...}, ...] list of dicts per gene
        self.__gene_info = [] # Stores info like {'name': 'gene1', 'type': 'DNA', 'length': 100, 'start': 1, 'end': 100}
        self.__all_taxons = []
        self.__concatenated_sequences = {} # {taxon_name: concatenated_sequence, ...}
        self.__partition_data = [] # [(gene_name, 'start-end', gene_type), ...]
        self.__statistics = {} # Dictionary holding various stats, including divergence


        self._parse_all_genes()
        self._collect_all_taxons()
        # _concatenate_sequences must run before _calculate_partition and _calculate_statistics
        # as it determines the final alignment length and gene positions.
        self.__concatenated_sequences = self._concatenate_sequences()
        # _calculate_partition requires gene_info which is updated by _concatenate_sequences
        self.__partition_data = self._calculate_partition() # Moved this definition above its call in __init__
        # _calculate_statistics requires concatenated_sequences and partition/gene_info
        self.__statistics = self._calculate_statistics()


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

    def _parse_all_genes(self) -> None:
        """
        Parses all raw gene file contents, detects format, and stores data and basic info.
        """
        for i, file_lines in enumerate(self.__raw_gene_contents):
            content_str = "".join(file_lines).strip()
            initial_gene_name = f"gene{i + 1}"
            parsed_data = {}
            determined_seq_type = 'Unknown'

            try:
                # Format detection
                is_nexus = content_str.lower().strip().startswith('#nexus') or 'begin data' in content_str.lower()
                is_genbank = content_str.lower().strip().startswith('locus') or 'origin' in content_str.lower() or 'version' in content_str.lower()
                is_fasta = content_str.strip().startswith('>')

                if is_nexus:
                    parsed_data = self._parse_nexus(file_lines)
                elif is_genbank:
                    parsed_data = self._parse_genbank(file_lines)
                elif is_fasta:
                    parsed_data = self._parse_fasta(file_lines)
                else:
                     # Try FASTA as a default if no specific format marker found
                    parsed_data = self._parse_fasta(file_lines)
                    if not parsed_data:
                         # If FASTA fails, maybe it's just plain sequence lines for one taxon?
                         # Simple attempt for a very basic format: assuming '>' not used, one taxon per file
                         if len(file_lines) > 0:
                             seq_content = "".join(line.strip() for line in file_lines if line.strip() and not line.strip().startswith('#') and not line.strip().startswith('[')).replace(' ', '').replace('\t','')
                             if seq_content:
                                  # Use a generic name as filename isn't available here easily
                                  dummy_taxon_name = f"Taxon_{initial_gene_name}" # Link taxon name to gene name if source is unknown format
                                  parsed_data = {dummy_taxon_name: seq_content}
                                  print(f"Info: Gene file {i+1} treated as simple single-sequence format.", file=sys.stderr)


                # Determine sequence type from parsed data if successful
                if parsed_data:
                    first_seq = None
                    for seq in parsed_data.values():
                         clean_seq = re.sub(r'[\s\-?]+', '', seq).upper()
                         if clean_seq:
                              first_seq = clean_seq
                              break

                    if first_seq:
                         dna_chars = set("ACGTU")
                         protein_chars = set("ACDEFGHIKLMNPQRSTVWY")
                         ambiguous_dna = set("RYSWKMBDHVN")
                         ambiguous_protein = set("BJZX")

                         all_dna_chars = dna_chars | ambiguous_dna | set('-?')
                         all_protein_chars = protein_chars | ambiguous_protein | set('-?')

                         seq_chars_present = set(first_seq)

                         is_potential_dna = seq_chars_present.issubset(all_dna_chars)
                         is_potential_protein = seq_chars_present.issubset(all_protein_chars)

                         if is_potential_dna and not (seq_chars_present & protein_chars):
                              determined_seq_type = 'DNA'
                         elif is_potential_protein and not (seq_chars_present & dna_chars):
                              determined_seq_type = 'Protein'
                         elif is_potential_dna and is_potential_protein:
                              if 'T' in seq_chars_present or 'U' in seq_chars_present:
                                   determined_seq_type = 'DNA'
                              elif seq_chars_present & set("FILPQEKRWYV"):
                                   determined_seq_type = 'Protein'
                              else:
                                   determined_seq_type = 'Unknown'
                         else:
                              determined_seq_type = 'Unknown'

                    if parsed_data:
                         self.__parsed_gene_data.append(parsed_data)
                         self.__gene_info.append({'name': initial_gene_name, 'type': determined_seq_type, 'length': 0, 'start': 0, 'end': 0})
                    else:
                         print(f"Warning: Gene file {i+1} parsed but yielded no sequences.", file=sys.stderr)


            except Exception as e:
                print(f"Error parsing gene file {i+1}: {e}", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)

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
                if match:
                    current_name = match.group(1).strip()
                    if not current_name: current_name = None
                else:
                     current_name = line[1:].strip()
                     if not current_name: current_name = None

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

        for raw_line in lines:
            line = raw_line.strip()
            if not line or line.startswith('#') or line.startswith('['):
                continue

            line_lower = line.lower()

            if 'begin data;' in line_lower:
                 sequences = {}
                 in_matrix = False
                 taxon_sequences_buffer = {}
                 current_taxon = None
                 continue

            if 'matrix;' in line_lower and not in_matrix:
                in_matrix = True
                taxon_sequences_buffer = {}
                current_taxon = None
                continue

            if 'end;' in line_lower:
                if in_matrix:
                     if current_taxon is not None and current_taxon in taxon_sequences_buffer:
                          sequences[current_taxon] = "".join(taxon_sequences_buffer[current_taxon])
                          del taxon_sequences_buffer[current_taxon]

                     for t, seq_parts in taxon_sequences_buffer.items():
                          if t not in sequences: # Handle interleaved format completion
                             sequences[t] = "".join(seq_parts)

                     in_matrix = False
                continue

            if in_matrix:
                match = re.match(r"^\s*(?:(['\"])(.*?)\1|(\S+))\s*(\S+)", line)
                if match:
                    quoted_name = match.group(2)
                    unquoted_name = match.group(3)
                    taxon_name = quoted_name if quoted_name is not None else unquoted_name
                    sequence_part = match.group(4).replace(';', '').strip()

                    if taxon_name:
                        if taxon_name not in taxon_sequences_buffer:
                            taxon_sequences_buffer[taxon_name] = []
                        taxon_sequences_buffer[taxon_name].append(sequence_part)

                        # If name changed and previous taxon had data, consider its sequence complete (sequential format)
                        if current_taxon is not None and current_taxon != taxon_name and current_taxon in taxon_sequences_buffer and taxon_sequences_buffer[current_taxon]:
                             if current_taxon not in sequences:
                                sequences[current_taxon] = "".join(taxon_sequences_buffer[current_taxon])
                             # Keep in buffer for interleaved check later
                        current_taxon = taxon_name

                elif current_taxon is not None:
                     sequence_part = line.replace(';', '').strip()
                     if sequence_part:
                          if current_taxon not in taxon_sequences_buffer:
                               taxon_sequences_buffer[current_taxon] = []
                          taxon_sequences_buffer[current_taxon].append(sequence_part)

        for t, seq_parts in taxon_sequences_buffer.items():
             if t not in sequences:
                  sequences[t] = "".join(seq_parts)

        return sequences

    def _parse_genbank(self, lines: list[str]) -> dict[str, str]:
        sequences = {}
        current_organism_name = None
        in_origin = False
        sequence_parts = []

        for line in lines:
            line_stripped = line.rstrip()
            if not line_stripped: continue

            if line_stripped == '//':
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

            if line_stripped.lower() == 'origin':
                in_origin = True
                sequence_parts = []
                continue

            if in_origin:
                sequence_parts.append(line_stripped)
                continue

        if current_organism_name is not None and sequence_parts:
             sequence = "".join(sequence_parts)
             clean_sequence = re.sub(r'[\s\d]+', '', sequence).upper()
             if clean_sequence:
                   sequences[current_organism_name] = clean_sequence

        return sequences

    def _collect_all_taxons(self) -> None:
        """
        Collects all unique taxon names from the parsed data.
        """
        all_taxons_set = set()
        for gene_data_dict in self.__parsed_gene_data:
            all_taxons_set.update(gene_data_dict.keys())
        self.__all_taxons = sorted(list(all_taxons_set))

    def _concatenate_sequences(self) -> dict[str, str]:
        """
        Concatenates sequences for each taxon, adding gaps for missing data.
        Updates gene_info with concatenated positions.
        """
        concatenated_sequences = {taxon: "" for taxon in self.__all_taxons}
        current_position = 1 # 1-based indexing

        # Iterate through the parsed data for each gene
        # We use the index 'i' to correspond to the gene info in __gene_info
        for i, gene_data_dict in enumerate(self.__parsed_gene_data):
            gene_length = 0
            if gene_data_dict:
                # Find the length of the first non-empty sequence for this gene
                first_non_empty_seq = next((seq for seq in gene_data_dict.values() if seq), "")
                gene_length = len(first_non_empty_seq)


            # Update the gene info with the determined length and position range
            # Check if the index 'i' is valid for __gene_info (should be if parsing was successful)
            if i < len(self.__gene_info):
                self.__gene_info[i]['length'] = gene_length
                self.__gene_info[i]['start'] = current_position
                self.__gene_info[i]['end'] = current_position + gene_length - 1
            else:
                 print(f"Warning: Mismatch between parsed data and gene info lists at index {i}.", file=sys.stderr)
                 # Create a dummy gene info entry if somehow missing (should ideally not happen)
                 self.__gene_info.append({'name': f'gene{i+1}_error', 'type': 'Unknown', 'length': gene_length, 'start': current_position, 'end': current_position + gene_length - 1})


            # Now, append the sequence (or gaps) for each taxon to the concatenated sequence
            for taxon in self.__all_taxons:
                sequence = gene_data_dict.get(taxon)

                if sequence is not None:
                    # If taxon has sequence data for this gene
                    if len(sequence) == gene_length:
                        # Sequence length matches expected length, append it
                        concatenated_sequences[taxon] += sequence
                    else:
                        # Sequence length is different from expected.
                        # Pad with gaps if shorter, or truncate if longer.
                        print(f"Warning: Taxon '{taxon}' sequence for gene '{self.__gene_info[i]['name']}' has unexpected length ({len(sequence)} vs {gene_length}). Padding/Truncating.", file=sys.stderr)
                        if len(sequence) < gene_length:
                            concatenated_sequences[taxon] += sequence + "-" * (gene_length - len(sequence))
                        else: # len(sequence) > gene_length
                            concatenated_sequences[taxon] += sequence[:gene_length]
                else:
                    # Taxon is missing from this gene, append gaps
                    concatenated_sequences[taxon] += "-" * gene_length

            # Update the current position for the next gene
            current_position += gene_length

        return concatenated_sequences

    def _calculate_partition(self) -> list[tuple[str, str, str]]:
        """
        Generates partition information based on processed gene_info.
        Positions are 1-based.
        This method must be defined *before* it's called in __init__.
        """
        partition = []
        for gene_info in self.__gene_info:
            gene_name = gene_info.get('name', 'UnknownGene')
            start = gene_info.get('start')
            end = gene_info.get('end')
            gene_type = gene_info.get('type', 'Unknown')

            # Only add partition if the gene contributed length to the alignment
            if start is not None and end is not None and start <= end:
                partition.append((gene_name, f"{start}-{end}", gene_type))
            elif start is not None and end is not None:
                 print(f"Warning: Invalid gene range for partition: {gene_name} ({start}-{end})", file=sys.stderr)
        return partition


    def _calculate_statistics(self) -> dict:
        """
        Calculates various statistics, including per-taxon per-gene divergence.
        """
        statistics = {}

        num_taxa = len(self.__all_taxons)
        statistics["Number of Taxa"] = num_taxa

        if num_taxa == 0 or not self.__concatenated_sequences:
             statistics["Total Length"] = 0
             statistics["Percentage Overall Missing Data (%)"] = 0.0
             statistics["Missing Data per Taxon (count)"] = {}
             statistics["Missing Data per Gene (count)"] = {}
             statistics["Gene Lengths (in concatenated alignment)"] = {}
             statistics["Number of Genes"] = len(self.__gene_info)
             statistics["Divergence Data"] = {}
             statistics["Taxon-Gene Matrix Sparsity (%)"] = 0.0
             return statistics

        total_length = len(next(iter(self.__concatenated_sequences.values()), "")) # Safely get length, defaults to 0 if no sequences
        statistics["Total Length"] = total_length

        total_cells = num_taxa * total_length
        total_missing_chars = 0
        missing_data_per_taxon = {}

        for taxon in self.__all_taxons:
            seq = self.__concatenated_sequences.get(taxon, "")
            missing_count = seq.count("-") + seq.count("?")
            missing_data_per_taxon[taxon] = missing_count
            total_missing_chars += missing_count

        statistics["Missing Data per Taxon (count)"] = missing_data_per_taxon

        missing_data_per_gene = {info['name']: 0 for info in self.__gene_info}
        for i, gene_info in enumerate(self.__gene_info):
             gene_name = gene_info['name']
             start = gene_info['start'] - 1 # Convert to 0-based index
             end = gene_info['end'] # End is inclusive in 1-based, so slice up to end (exclusive)

             if gene_name in missing_data_per_gene: # Ensure gene_name is expected
                  if start >= 0 and end >= start and total_length >= end:
                       for taxon in self.__all_taxons:
                            seq = self.__concatenated_sequences.get(taxon, "")
                            if seq and len(seq) > start:
                                 seq_segment = seq[start:min(end, len(seq))]
                                 missing_count_segment = seq_segment.count("-") + seq_segment.count("?")
                                 missing_data_per_gene[gene_name] += missing_count_segment
                  elif gene_info.get('length', 0) > 0:
                       print(f"Warning: Invalid gene range calculated for '{gene_name}' ({gene_info['start']}-{gene_info['end']}). Skipping missing data count for this gene.", file=sys.stderr)


        statistics["Missing Data per Gene (count)"] = missing_data_per_gene

        if total_cells > 0:
            percentage_missing_overall = (total_missing_chars / total_cells) * 100
            statistics["Percentage Overall Missing Data (%)"] = round(percentage_missing_overall, 2)
        else:
            statistics["Percentage Overall Missing Data (%)"] = 0.0

        num_expected_gene_segments = len(self.__all_taxons) * len(self.__gene_info)
        num_fully_missing_gene_segments = 0
        for taxon in self.__all_taxons:
             seq = self.__concatenated_sequences.get(taxon, "")
             for gene_info in self.__gene_info:
                  start_0based = gene_info['start'] - 1
                  end_0based = gene_info['end']
                  gene_segment_length = gene_info.get('length', 0)

                  if start_0based >= 0 and end_0based >= start_0based and seq and len(seq) >= end_0based:
                       gene_segment = seq[start_0based:end_0based]
                       if gene_segment_length > 0 and len(gene_segment) == gene_segment_length and all(char in ('-', '?') for char in gene_segment):
                            num_fully_missing_gene_segments += 1

        if num_expected_gene_segments > 0:
             sparsity_percentage = (num_fully_missing_gene_segments / num_expected_gene_segments) * 100
             statistics["Taxon-Gene Matrix Sparsity (%)"] = round(sparsity_percentage, 2)
        else:
             statistics["Taxon-Gene Matrix Sparsity (%)"] = 0.0


        gene_lengths_in_concat = {info['name']: info.get('length', 0) for info in self.__gene_info}
        statistics["Gene Lengths (in concatenated alignment)"] = gene_lengths_in_concat
        statistics["Number of Genes"] = len(self.__gene_info)


        divergence_data = self._calculate_per_taxon_per_gene_stats(self.__concatenated_sequences, list(self.__gene_info), list(self.__all_taxons))
        statistics["Divergence Data"] = divergence_data


        return statistics

    def _calculate_per_taxon_per_gene_stats(self, concatenated_sequences: dict[str, str], gene_info: list[dict], all_taxons: list[str]) -> dict:
        """
        Calculates difference statistics for each taxon relative to the first taxon
        (assumed reference) for each gene segment and overall.
        """
        divergence_data = {}

        if not all_taxons or len(all_taxons) < 1 or not concatenated_sequences or not gene_info:
            for taxon in all_taxons:
                 divergence_data[taxon] = {
                     'Total score': 0.0,
                     'No of charsets': 0,
                 }
            return divergence_data


        reference_taxon = all_taxons[0]
        ref_seq = concatenated_sequences.get(reference_taxon, "")
        total_alignment_length = len(ref_seq)

        # Initialize data structure for each taxon, including N/A for all gene scores initially
        for taxon in all_taxons:
             divergence_data[taxon] = {
                  'Total score': 0.0,
                  'No of charsets': 0,
             }
             for gene in gene_info:
                 divergence_data[taxon][gene['name']] = "N/A" # Default for gene scores


        # Calculate per-gene differences and count non-gap charsets
        for gene in gene_info:
             gene_name = gene['name']
             start_0based = gene['start'] - 1
             end_0based = gene['end']

             if start_0based < 0 or end_0based <= start_0based or end_0based > total_alignment_length:
                  print(f"Warning: Invalid gene range for '{gene_name}' ({gene['start']}-{gene['end']}) within alignment length {total_alignment_length}. Skipping divergence calculation for this gene.", file=sys.stderr)
                  continue

             ref_segment = ref_seq[start_0based:end_0based]
             gene_segment_length = len(ref_segment)

             if gene_segment_length == 0:
                  print(f"Warning: Gene '{gene_name}' has calculated length 0. Skipping divergence calculation for this gene.", file=sys.stderr)
                  continue

             for taxon in all_taxons:
                  if taxon == reference_taxon:
                       divergence_data[taxon][gene_name] = f"0% #0 (0.00%)"
                       divergence_data[taxon]['No of charsets'] = len([g for g in gene_info if g.get('length', 0) > 0])
                       continue


                  taxon_seq = concatenated_sequences.get(taxon, "")
                  if len(taxon_seq) < end_0based:
                       # This shouldn't happen if concatenate works correctly, but defensive slice
                       taxon_segment = taxon_seq[start_0based:] + "-" * (end_0based - len(taxon_seq))
                  else:
                       taxon_segment = taxon_seq[start_0based:end_0based]


                  diff_count = 0
                  comparable_sites = 0
                  taxon_has_non_gap_in_gene_segment = False

                  for j in range(gene_segment_length):
                       ref_char = ref_segment[j]
                       taxon_char = taxon_segment[j]

                       if taxon_char not in ('-', '?'):
                            taxon_has_non_gap_in_gene_segment = True

                       if ref_char not in ('-', '?') and taxon_char not in ('-', '?'):
                            comparable_sites += 1
                            if ref_char.upper() != taxon_char.upper():
                                 diff_count += 1

                  # Calculate percentage based on gene segment length, matching the image format logic
                  percentage = (diff_count / gene_segment_length) * 100 if gene_segment_length > 0 else 0.0

                  leading_percent = int(round(percentage))
                  if leading_percent > 99: leading_percent = 99
                  if leading_percent < 0: leading_percent = 0

                  divergence_data[taxon][gene_name] = f"{leading_percent}% #{diff_count} ({percentage:.2f}%)"

                  # Increment No of charsets if taxon has any non-gap data in this specific gene segment
                  if taxon_has_non_gap_in_gene_segment:
                       divergence_data[taxon]['No of charsets'] += 1


        # Calculate Total Score per taxon (Average of gene percentages where taxon had non-gap data)
        for taxon in all_taxons:
             if taxon == reference_taxon:
                  divergence_data[taxon]['Total score'] = 0.0
                  continue

             total_percentage_sum = 0.0
             genes_counted_for_total_score = 0

             for gene in gene_info:
                  gene_name = gene['name']
                  gene_score_str = divergence_data[taxon].get(gene_name)

                  if gene_score_str and "N/A" not in gene_score_str:
                       match = re.search(r'\((\d+\.?\d*)%\)', gene_score_str) # Handle cases like (0%) or (10.5%)
                       if match:
                            try:
                                 percentage_value = float(match.group(1))
                                 total_percentage_sum += percentage_value
                                 genes_counted_for_total_score += 1
                            except ValueError:
                                 print(f"Warning: Could not parse percentage from gene score string for taxon '{taxon}', gene '{gene_name}': {gene_score_str}", file=sys.stderr)
                       elif gene_score_str.startswith("0%") and re.search(r'\(\d+\.?\d*%\)', gene_score_str): # Handle exact "0% #X (Y.YY%)" case
                             # Percentage is 0.0
                             genes_counted_for_total_score += 1 # Count this gene even if % is 0


             average_percentage = (total_percentage_sum / genes_counted_for_total_score) if genes_counted_for_total_score > 0 else 0.0
             divergence_data[taxon]['Total score'] = round(average_percentage, 2)


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
    ]
    dummy_fasta_gene2 = [
        ">TaxonA\nTTTT----GGGG", # Length 12
        ">TaxonB\nCCCCAAAA----",
        ">TaxonD\nGGGGGGGGGGGG"
    ]
    dummy_fasta_gene3 = [
        ">TaxonA\nAAAAAAAAAA", # Length 10
        ">TaxonB\nTTTTTTTTTT",
        ">TaxonC\nCCCCCCCCCC"
    ]


    raw_contents = [
        dummy_fasta_gene1,
        dummy_fasta_gene2,
        dummy_fasta_gene3,
    ]

    concatenator = SequenceConcatenator(raw_contents)

    concatenated = concatenator.get_concatenated_sequences()
    stats = concatenator.get_statistics()
    partition = concatenator.get_partition()

    print("\nConcatenated Sequences:")
    for taxon, seq in concatenated.items():
        print(f">{taxon}\n{seq}")

    print("\nStatistics:")
    import json
    class SetEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, set):
                return list(obj)
            return json.JSONEncoder.default(self, obj)

    print(json.dumps(stats, indent=2, cls=SetEncoder))

    print("\nPartition:")
    for gene_name, pos_range, gene_type in partition:
        print(f"Gene: {gene_name}, Range: {pos_range}, Type: {gene_type}")