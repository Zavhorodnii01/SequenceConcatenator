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
        self.__partition_data = self._calculate_partition()
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
        Uses the gene name provided by the frontend.
        """
        self.__parsed_gene_data = [] # Clear existing data on re-parse
        self.__gene_info = []      # Clear existing info on re-parse

        for i, file_lines in enumerate(self.__raw_gene_contents):
            content_str = "".join(file_lines).strip()
            # Get the gene name from the list passed from the frontend
            # Use a fallback if somehow out of index (should not happen if init check passes)
            current_gene_name = self.__frontend_gene_names[i] if i < len(self.__frontend_gene_names) else f"gene_index_{i+1}"

            parsed_data = {}
            determined_seq_type = 'Unknown'

            try:
                # Format detection (case-insensitive)
                is_nexus = content_str.lower().strip().startswith('#nexus') or 'begin data' in content_str.lower()
                is_genbank = content_str.lower().strip().startswith('locus') or 'origin' in content_str.lower() or 'version' in content_str.lower() or 'ncbi' in content_str.lower()
                is_fasta = content_str.strip().startswith('>') # FASTA header must start with '>'

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
                         # If FASTA fails and it's not explicitly NEXUS/GenBank,
                         # maybe it's just plain sequence lines for one taxon?
                         # This is a highly speculative format, handle with care.
                         # Only proceed if there's actual content.
                         non_empty_lines = [line.strip() for line in file_lines if line.strip() and not line.strip().startswith('#') and not line.strip().startswith('[')]
                         if non_empty_lines:
                              seq_content = "".join(non_empty_lines).replace(' ', '').replace('\t','')
                              if seq_content:
                                  # Use the gene name provided by the frontend as a basis for a dummy taxon name
                                  # This handles the edge case where a file only contains sequence data without a header
                                  dummy_taxon_name = f"Taxon_for_{current_gene_name}_File"
                                  print(f"Info: Gene file index {i+1} ('{current_gene_name}') contains data but no recognized header. Treating as simple single-sequence for taxon '{dummy_taxon_name}'.", file=sys.stderr)
                                  parsed_data = {dummy_taxon_name: seq_content}
                                  if not seq_content:
                                      print(f"Warning: Gene file index {i+1} ('{current_gene_name}') had lines but yielded no sequence data after cleaning.", file=sys.stderr)
                                      parsed_data = {} # Treat as empty if no sequence resulted


                # Determine sequence type from parsed data if successful
                # This logic block remains mostly the same, operating on parsed_data
                if parsed_data:
                    first_seq = None
                    for seq in parsed_data.values():
                         clean_seq = re.sub(r'[\s\-?]+', '', seq).upper()
                         if clean_seq:
                              first_seq = clean_seq
                              break # Found the first non-empty sequence

                    if first_seq:
                         dna_chars = set("ACGTU")
                         protein_chars = set("ACDEFGHIKLMNPQRSTVWY")
                         ambiguous_dna = set("RYSWKMBDHVN")
                         ambiguous_protein = set("BJZX") # Includes B (Asp/Asn), J (Leu/Ile), Z (Glu/Gln), X (Any)

                         # Include gap/missing characters in the allowed set for parsing,
                         # but not for determining the *primary* type unless it's all gaps.
                         all_dna_chars_strict = dna_chars | ambiguous_dna
                         all_protein_chars_strict = protein_chars | ambiguous_protein


                         seq_chars_present = set(first_seq) # Chars in the *first non-empty sequence*

                         is_potential_dna = seq_chars_present.issubset(all_dna_chars_strict)
                         is_potential_protein = seq_chars_present.issubset(all_protein_chars_strict)

                         # Refined type determination logic
                         if not seq_chars_present:
                             determined_seq_type = 'Unknown' # Empty sequence
                         elif is_potential_dna and not (seq_chars_present & protein_chars):
                              # Looks like pure DNA characters (including ambiguous DNA, excluding standard protein)
                              determined_seq_type = 'DNA'
                         elif is_potential_protein and not (seq_chars_present & dna_chars):
                              # Looks like pure Protein characters (including ambiguous protein, excluding standard DNA)
                              determined_seq_type = 'Protein'
                         elif is_potential_dna and is_potential_protein:
                              # Contains characters common to both (A, C, G, T -> Gly, Ala, Asp, Cys etc.) or N/X/?.
                              # Look for discriminating characters.
                              if seq_chars_present & set("TU"): # T (DNA), U (RNA) are strong indicators of nucleic acid
                                   determined_seq_type = 'DNA'
                              elif seq_chars_present & set("FILPQEKRWYV"): # Amino acids less common or not in DNA alphabet
                                   determined_seq_type = 'Protein'
                              else:
                                   # Still ambiguous (e.g., only A, C, G, -, ?). Assume DNA as it's more common for concatenation? Or Unknown.
                                   # Stick with Unknown if ambiguous.
                                   determined_seq_type = 'Unknown'
                         else:
                              determined_seq_type = 'Unknown' # Contains characters not in either set

                    if parsed_data:
                         self.__parsed_gene_data.append(parsed_data)
                         # Use the current_gene_name obtained from the frontend list
                         self.__gene_info.append({'name': current_gene_name, 'type': determined_seq_type, 'length': 0, 'start': 0, 'end': 0})
                    else:
                         print(f"Warning: Gene file index {i+1} ('{current_gene_name}') parsed but yielded no sequences.", file=sys.stderr)


            except Exception as e:
                print(f"Error parsing gene file index {i+1} ('{current_gene_name}'): {e}", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)

        # After parsing all genes, filter out genes that resulted in no parsed data
        # This ensures gene_info and parsed_gene_data stay aligned and only contain valid genes
        valid_genes_data = []
        valid_gene_info = []
        for i, gene_data in enumerate(self.__parsed_gene_data):
             # Check if gene_data is not empty and has sequences
             if gene_data and any(seq.strip() for seq in gene_data.values()):
                  valid_genes_data.append(gene_data)
                  # Find the corresponding info from the original __gene_info list
                  # This lookup assumes __gene_info still matches the original __parsed_gene_data list indices before filtering
                  if i < len(self.__gene_info):
                       valid_gene_info.append(self.__gene_info[i])
                  else:
                       # Fallback if index is somehow off (should not happen)
                       print(f"Internal Warning: Mismatch indexing gene info after filtering.", file=sys.stderr)


        self.__parsed_gene_data = valid_genes_data
        self.__gene_info = valid_gene_info
        # Update frontend gene names list to match filtered results?
        # No, the frontend list tracks what the user ADDED. The backend works with what it SUCCESSFULLY PARSED.
        # The frontend should probably not rely on the backend's internal __gene_info list for its gene listbox.
        # The frontend's gene_names list and raw_gene_contents list are the source of truth for the input tab.
        # The backend processes these and produces results based on successful parsing.
        # The mismatch in counts should be handled gracefully by the UI (e.g., statistics showing fewer genes processed).
        # The current logic in on_submit already catches if _concatenated_sequences is empty, which covers parsing failures.


    def _parse_fasta(self, lines: list[str]) -> dict[str, str]:
        sequences = {}
        current_name = None
        current_seq_lines = []

        # Process lines, ignoring empty lines and comments
        processed_lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]

        for line in processed_lines:
            if line.startswith(">"):
                # Save the previous sequence if any
                if current_name is not None and current_seq_lines:
                    sequences[current_name] = "".join(current_seq_lines)

                # Parse the new header
                match = re.match(r'>\s*(.+)', line)
                if match:
                    current_name = match.group(1).strip()
                    if not current_name:
                         # Handle case of '>' with no name
                         current_name = f"UnnamedTaxon_{len(sequences) + 1}"
                         print(f"Warning: Found FASTA header '>' with no name. Assigning default name '{current_name}'.", file=sys.stderr)
                else:
                     # Should not happen with startswith('>'), but as a fallback
                     current_name = line[1:].strip()
                     if not current_name:
                          current_name = f"UnnamedTaxon_{len(sequences) + 1}"
                          print(f"Warning: Found FASTA header line '{line}' with no name. Assigning default name '{current_name}'.", file=sys.stderr)

                current_seq_lines = [] # Reset sequence lines for the new taxon
            elif current_name is not None:
                # Append sequence lines to the current taxon, clean out whitespace and digits
                clean_seq_line = re.sub(r'[\s\d]+', '', line)
                current_seq_lines.append(clean_seq_line)

        # Save the last sequence after the loop finishes
        if current_name is not None and current_seq_lines:
            sequences[current_name] = "".join(current_seq_lines)

        # Clean sequences: convert to uppercase and remove any remaining illegal characters if necessary (optional)
        # Let's keep the cleaning minimal to just whitespace/digits during parsing
        # Type detection happens later based on characters present

        return sequences

    def _parse_nexus(self, lines: list[str]) -> dict[str, str]:
        sequences = {}
        in_matrix = False
        taxon_sequences_buffer = {}
        current_taxon = None
        inside_block_level = 0 # Track nested BEGIN/END blocks

        for raw_line in lines:
            line = raw_line.strip()
            if not line: # Ignore empty lines
                continue

            # Simple comment handling (NEXUS uses [])
            if '[' in line and ']' in line:
                 line = re.sub(r'\[.*?\]', '', line).strip()
                 if not line: continue # Line was only a comment

            line_lower = line.lower()

            # Handle BEGIN/END blocks to only parse MATRIX inside DATA
            if line_lower.startswith('begin'):
                 inside_block_level += 1
                 if ' data;' in line_lower and inside_block_level == 1:
                      # Reset data parsing state when entering DATA block
                      sequences = {}
                      in_matrix = False
                      taxon_sequences_buffer = {}
                      current_taxon = None
                 continue

            if line_lower.startswith('end;'):
                 if in_matrix:
                      # If we were inside the matrix block, process the last buffered taxon sequence
                      if current_taxon is not None and current_taxon in taxon_sequences_buffer:
                           sequences[current_taxon] = "".join(taxon_sequences_buffer[current_taxon])
                           del taxon_sequences_buffer[current_taxon] # Clear buffer for this taxon

                     # After leaving MATRIX, add any remaining buffered sequences (handles interleaved)
                      for t, seq_parts in list(taxon_sequences_buffer.items()): # Iterate over a copy
                           if t not in sequences:
                              sequences[t] = "".join(seq_parts)
                           del taxon_sequences_buffer[t] # Clear buffer

                      in_matrix = False

                 inside_block_level = max(0, inside_block_level - 1) # Decrement block level
                 continue # Skip to next line after processing END

            # Look for MATRIX inside the DATA block
            if inside_block_level == 1 and 'matrix;' in line_lower:
                in_matrix = True
                taxon_sequences_buffer = {} # Reset buffer for new matrix block
                current_taxon = None
                continue

            if in_matrix:
                # Parse lines within the matrix block
                # Use a regex that handles quoted or unquoted names followed by sequence
                match = re.match(r"^\s*(?:(['\"])(.*?)\1|([^'\s]+))\s*(\S+)\s*;", line)
                if not match:
                     # Try again without trailing semicolon, some files omit it on every line
                     match = re.match(r"^\s*(?:(['\"])(.*?)\1|([^'\s]+))\s*(\S+)", line)

                if match:
                    quoted_name = match.group(2)
                    unquoted_name = match.group(3)
                    taxon_name = quoted_name if quoted_name is not None else unquoted_name
                    sequence_part = match.group(4).strip() # Ensure no trailing/leading whitespace

                    if taxon_name:
                        # If we encounter a new taxon name, save the sequence of the previous one
                        if current_taxon is not None and current_taxon != taxon_name and current_taxon in taxon_sequences_buffer:
                             # If the buffer for the previous taxon has parts, this indicates a new sequential block or interleaved block line
                             # Combine buffered parts and add to sequences dictionary if not already there
                             if taxon_sequences_buffer[current_taxon]:
                                  if current_taxon not in sequences:
                                       sequences[current_taxon] = "".join(taxon_sequences_buffer[current_taxon])
                                  # Keep in buffer for interleaved check until END; DATA is hit? Or clear?
                                  # Clearing seems safer for sequential. For interleaved, need to keep.
                                  # Let's keep in buffer until END; DATA.

                        current_taxon = taxon_name # Update the current taxon being processed

                        # Append the sequence part to the buffer for this taxon
                        if taxon_name not in taxon_sequences_buffer:
                            taxon_sequences_buffer[taxon_name] = []
                        taxon_sequences_buffer[taxon_name].append(sequence_part)

                elif current_taxon is not None:
                     # Handle lines that might just be sequence data in interleaved format
                     sequence_part = line.replace(';', '').strip()
                     if sequence_part:
                          if current_taxon not in taxon_sequences_buffer:
                               # This case implies a line of sequence without a preceding taxon name,
                               # but we assume it belongs to the *current_taxon*. This is the nature of interleaved.
                               # If current_taxon is None here, it's likely a parsing error or malformed file.
                               pass # Already handled by current_taxon check
                          else:
                               taxon_sequences_buffer[current_taxon].append(sequence_part)
                     # If line is empty or just ';', ignore within matrix unless it's the last line?

        # After the loop, add any sequences remaining in the buffer (important for interleaved)
        for t, seq_parts in taxon_sequences_buffer.items():
             if t not in sequences:
                  sequences[t] = "".join(seq_parts)

        # Clean sequences: remove any remaining whitespace, newlines, or digits from concatenated parts
        cleaned_sequences = {name: re.sub(r'[\s\d]+', '', seq).upper() for name, seq in sequences.items() if seq.strip()}
        return cleaned_sequences


    def _parse_genbank(self, lines: list[str]) -> dict[str, str]:
        sequences = {}
        current_organism_name = None
        in_origin = False
        sequence_parts = []

        for line in lines:
            line_stripped = line.rstrip() # Keep trailing whitespace for keywords
            if not line_stripped.strip(): continue # Ignore empty lines

            # End of record marker
            if line_stripped.strip() == '//':
                if current_organism_name is not None and sequence_parts:
                    sequence = "".join(sequence_parts)
                    clean_sequence = re.sub(r'[\s\d]+', '', sequence).upper() # Remove whitespace and digits from sequence

                    if clean_sequence: # Only add if sequence is not empty after cleaning
                         sequences[current_organism_name] = clean_sequence
                         # print(f"Debug: Parsed GenBank sequence for '{current_organism_name}' (length {len(clean_sequence)})") # Debugging line

                # Reset state for the next record
                current_organism_name = None
                in_origin = False
                sequence_parts = []
                continue

            # Organism line (usually before ORIGIN)
            org_match = re.match(r'^\s{2,}ORGANISM\s+(.*?)\s*\.?$', line_stripped, re.IGNORECASE)
            if org_match:
                current_organism_name = org_match.group(1).strip()
                in_origin = False # Ensure we are not in ORIGIN block when finding organism
                sequence_parts = [] # Reset parts for the new record
                # print(f"Debug: Found Organism: '{current_organism_name}'") # Debugging line
                continue

            # ORIGIN line marks the start of sequence data
            if line_stripped.strip().lower() == 'origin':
                in_origin = True
                sequence_parts = [] # Reset parts for the sequence
                # print("Debug: Entered ORIGIN block.") # Debugging line
                continue

            # Lines within the ORIGIN block contain sequence data
            if in_origin:
                # Append the raw line within ORIGIN. Cleaning happens later.
                sequence_parts.append(line_stripped)
                # print(f"Debug: Added line to sequence_parts: '{line_stripped}'") # Debugging line
                continue

            # Additional check: Sometimes the taxon name might be in a FEATURES block, e.g., /organism="..."
            # This is less standard for the primary taxon name, but can be a fallback.
            # Let's add a simple check, prioritizing the main ORGANISM line.
            if current_organism_name is None:
                 org_feature_match = re.search(r'/organism="([^"]+)"', line_stripped, re.IGNORECASE)
                 if org_feature_match:
                      current_organism_name = org_feature_match.group(1).strip()
                      # print(f"Debug: Found organism in feature: '{current_organism_name}'") # Debugging line
                      # Note: This doesn't reset sequence_parts etc., as it might be mid-file.
                      # It just sets the name for the *current* record if not already set.


        # After the loop, process the last record if the file didn't end with //
        # (This handles files that are just a single GenBank record without a final //)
        if current_organism_name is not None and sequence_parts:
             sequence = "".join(sequence_parts)
             clean_sequence = re.sub(r'[\s\d]+', '', sequence).upper() # Remove whitespace and digits

             if clean_sequence:
                   sequences[current_organism_name] = clean_sequence
                   # print(f"Debug: Processed last GenBank sequence for '{current_organism_name}' (length {len(clean_sequence)})") # Debugging line


        # Filter out entries with empty sequences after cleaning
        cleaned_sequences = {name: seq for name, seq in sequences.items() if seq}
        return cleaned_sequences


    def _collect_all_taxons(self) -> None:
        """
        Collects all unique taxon names from the parsed data.
        """
        all_taxons_set = set()
        for gene_data_dict in self.__parsed_gene_data:
            all_taxons_set.update(gene_data_dict.keys())
        self.__all_taxons = sorted(list(all_taxons_set))
        # print(f"Debug: Collected {len(self.__all_taxons)} unique taxons: {self.__all_taxons}") # Debugging


    def _concatenate_sequences(self) -> dict[str, str]:
        """
        Concatenates sequences for each taxon across all parsed genes.
        Adds gaps for missing data or sequences of unexpected length.
        Updates gene_info with concatenated positions.
        """
        concatenated_sequences = {taxon: "" for taxon in self.__all_taxons}
        current_position = 1 # 1-based indexing for partition/gene_info

        # Ensure gene_info is cleared and rebuilt based on currently parsed data
        # This should already happen in _parse_all_genes, but defensive clear
        # self.__gene_info = [] # No, this should be populated by _parse_all_genes

        temp_gene_info = [] # Build the new gene_info list based on actual sequence data length

        # Iterate through the parsed data for each gene
        # We use the index 'i' to correspond to the parsed gene data dictionary
        for i, gene_data_dict in enumerate(self.__parsed_gene_data):
            # Get the corresponding gene name from the filtered __gene_info list
            # This assumes __gene_info and __parsed_gene_data are aligned after filtering in _parse_all_genes
            if i >= len(self.__gene_info):
                 print(f"Internal Error: Mismatch between __parsed_gene_data and __gene_info lists at index {i} during concatenation.", file=sys.stderr)
                 # Skip this gene if info is missing, or use a fallback name
                 gene_name = f"Gene_Error_{i+1}"
                 gene_type = "Unknown"
                 gene_length = 0
                 # Attempt to find the first non-empty sequence length anyway
                 if gene_data_dict:
                     first_non_empty_seq = next((seq for seq in gene_data_dict.values() if seq), "")
                     gene_length = len(first_non_empty_seq)

            else:
                 gene_name = self.__gene_info[i].get('name', f'UnnamedGene_{i+1}')
                 gene_type = self.__gene_info[i].get('type', 'Unknown')
                 # Determine the *actual* length of this gene segment from the parsed data
                 # This is the length all sequences for this gene should ideally have after parsing/cleaning
                 gene_length = 0
                 if gene_data_dict:
                     # Find the length of the first non-empty sequence for this gene
                     first_non_empty_seq = next((seq for seq in gene_data_dict.values() if seq), "")
                     gene_length = len(first_non_empty_seq)
                     # Use the max length across all taxons for this gene? No, usually consensus length is best.
                     # The current logic uses the first non-empty sequence length. Let's stick to that.


            if gene_length > 0:
                 # Update the gene info with the determined length and position range
                 gene_start = current_position
                 gene_end = current_position + gene_length - 1
                 temp_gene_info.append({'name': gene_name, 'type': gene_type, 'length': gene_length, 'start': gene_start, 'end': gene_end})

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
                 # If gene_length is 0 (e.g., no sequences parsed for this gene, or all were empty),
                 # this gene contributes 0 length to the alignment, so current_position doesn't change.
                 # We still add a gene_info entry with length 0, although it might be skipped for partition later.
                 print(f"Info: Gene '{gene_name}' has effective length 0 after parsing. Not added to concatenated alignment.", file=sys.stderr)
                 temp_gene_info.append({'name': gene_name, 'type': gene_type, 'length': 0, 'start': current_position, 'end': current_position -1}) # End before start


        # Update the official gene_info list with the calculated lengths and positions
        self.__gene_info = temp_gene_info

        # Final check: ensure all concatenated sequences have the same total length
        if self.__all_taxons:
            expected_total_length = sum(info.get('length', 0) for info in self.__gene_info)
            for taxon, seq in concatenated_sequences.items():
                 if len(seq) != expected_total_length:
                      print(f"Internal Error: Taxon '{taxon}' concatenated sequence has unexpected total length ({len(seq)} vs {expected_total_length}). This indicates a bug in concatenation logic.", file=sys.stderr)
                      # This could lead to errors later, consider raising an exception or trying to fix the length
                      # For now, just print error, the UI might handle inconsistent sequence lengths poorly.

        return concatenated_sequences

    def _calculate_partition(self) -> list[tuple[str, str, str]]:
        """
        Generates partition information based on processed gene_info.
        Positions are 1-based.
        This method must be defined *before* it's called in __init__.
        """
        partition = []
        # Use the gene_info list updated by _concatenate_sequences
        for gene_info in self.__gene_info:
            gene_name = gene_info.get('name', 'UnknownGene')
            start = gene_info.get('start')
            end = gene_info.get('end')
            gene_type = gene_info.get('type', 'Unknown')

            # Only add partition if the gene contributed length to the alignment
            if start is not None and end is not None and start <= end:
                partition.append((gene_name, f"{start}-{end}", gene_type))
            # No warning here for length 0 genes, they are intentionally skipped from partition block

        return partition


    def _calculate_statistics(self) -> dict:
        """
        Calculates various statistics, including per-taxon per-gene divergence.
        This method uses the potentially filtered self.__gene_info list.
        """
        statistics = {}

        # Calculate statistics based on the FINAL concatenated sequences and gene info
        num_taxa = len(self.__concatenated_sequences) # Use count from concatenated sequences (reflects common taxa)
        statistics["Number of Taxa"] = num_taxa
        statistics["Number of Genes"] = len([g for g in self.__gene_info if g.get('length', 0) > 0]) # Count genes that contributed sequence

        # Get the actual list of taxons present in the concatenated data
        taxons_in_concat = sorted(list(self.__concatenated_sequences.keys()))

        if num_taxa == 0 or not self.__concatenated_sequences:
             statistics["Total Length"] = 0
             statistics["Percentage Overall Missing Data (%)"] = 0.0
             statistics["Missing Data per Taxon (count)"] = {}
             statistics["Missing Data per Gene (count)"] = {}
             statistics["Gene Lengths (in concatenated alignment)"] = {}
             statistics["Taxon-Gene Matrix Sparsity (%)"] = 0.0
             statistics["Divergence Data"] = {}
             return statistics

        total_length = len(next(iter(self.__concatenated_sequences.values()), "")) # Safely get length
        statistics["Total Length"] = total_length

        total_cells = num_taxa * total_length
        total_missing_chars = 0
        missing_data_per_taxon = {}

        for taxon in taxons_in_concat:
            seq = self.__concatenated_sequences.get(taxon, "")
            missing_count = seq.count("-") + seq.count("?")
            missing_data_per_taxon[taxon] = missing_count
            total_missing_chars += missing_count

        statistics["Missing Data per Taxon (count)"] = missing_data_per_taxon


        # Calculate missing data per gene based on concatenated sequences
        missing_data_per_gene = {info['name']: 0 for info in self.__gene_info if info.get('length', 0) > 0} # Only include genes with length > 0
        for gene_info in self.__gene_info:
             gene_name = gene_info['name']
             start_0based = gene_info['start'] - 1 # Convert to 0-based index
             end_0based = gene_info['end'] # End is inclusive in 1-based, so slice up to end (exclusive)
             gene_segment_length = gene_info.get('length', 0)

             if gene_segment_length > 0 and gene_name in missing_data_per_gene: # Ensure gene was included and has length
                  if start_0based >= 0 and end_0based >= start_0based and total_length >= end_0based:
                       for taxon in taxons_in_concat:
                            seq = self.__concatenated_sequences.get(taxon, "")
                            if seq and len(seq) >= end_0based: # Ensure seq is long enough
                                 seq_segment = seq[start_0based:end_0based]
                                 missing_count_segment = seq_segment.count("-") + seq_segment.count("?")
                                 missing_data_per_gene[gene_name] += missing_count_segment
                            # else: # Taxon seq is shorter than end_0based - should be handled by concatenation padding, but defensive check
                            #     print(f"Internal Warning: Taxon sequence length unexpected for '{taxon}' gene '{gene_name}'.")

                  # else: # Invalid range warning is already given in concatenation
                  #     print(f"Warning: Invalid gene range calculated for '{gene_name}' ({gene_info['start']}-{gene_info['end']}). Skipping missing data count for this gene.", file=sys.stderr)


        statistics["Missing Data per Gene (count)"] = missing_data_per_gene

        if total_cells > 0:
            percentage_missing_overall = (total_missing_chars / total_cells) * 100
            statistics["Percentage Overall Missing Data (%)"] = round(percentage_missing_overall, 2)
        else:
            statistics["Percentage Overall Missing Data (%)"] = 0.0

        # Calculate Taxon-Gene Matrix Sparsity
        num_genes_with_length = len([g for g in self.__gene_info if g.get('length', 0) > 0])
        num_expected_gene_segments = num_taxa * num_genes_with_length
        num_fully_missing_gene_segments = 0

        if num_genes_with_length > 0 and num_taxa > 0:
             for taxon in taxons_in_concat:
                  seq = self.__concatenated_sequences.get(taxon, "")
                  for gene_info in self.__gene_info:
                       gene_segment_length = gene_info.get('length', 0)
                       if gene_segment_length > 0: # Only count segments for genes that contributed length
                            start_0based = gene_info['start'] - 1
                            end_0based = gene_info['end']

                            if seq and len(seq) >= end_0based:
                                 gene_segment = seq[start_0based:end_0based]
                                 if len(gene_segment) == gene_segment_length and all(char in ('-', '?') for char in gene_segment):
                                      num_fully_missing_gene_segments += 1
                            else:
                                 # This case indicates the sequence for this taxon is unexpectedly short,
                                 # but the concatenation logic *should* have padded it. If it didn't, it's a bug.
                                 # Assuming padding works, a sequence exists and is >= end_0based.
                                 # If it reaches here, it implies a logical error or a taxon completely missing from concat data (already filtered by taxons_in_concat).
                                 pass # Should not happen if concatenation is correct


        if num_expected_gene_segments > 0:
             sparsity_percentage = (num_fully_missing_gene_segments / num_expected_gene_segments) * 100
             statistics["Taxon-Gene Matrix Sparsity (%)"] = round(sparsity_percentage, 2)
        else:
             statistics["Taxon-Gene Matrix Sparsity (%)"] = 0.0


        gene_lengths_in_concat = {info['name']: info.get('length', 0) for info in self.__gene_info}
        statistics["Gene Lengths (in concatenated alignment)"] = gene_lengths_in_concat


        # Calculate divergence data using the final list of taxons and gene info
        divergence_data = self._calculate_per_taxon_per_gene_stats(self.__concatenated_sequences, self.__gene_info, taxons_in_concat)
        statistics["Divergence Data"] = divergence_data


        return statistics

    def _calculate_per_taxon_per_gene_stats(self, concatenated_sequences: dict[str, str], gene_info: list[dict], taxons_in_concat: list[str]) -> dict:
        """
        Calculates difference statistics for each taxon relative to the first taxon
        (assumed reference) for each gene segment and overall.
        Uses actual gene names from gene_info.
        """
        divergence_data = {}

        # Use the list of taxons actually present in the concatenated alignment
        taxons = taxons_in_concat

        if not taxons or len(taxons) < 1 or not concatenated_sequences or not gene_info:
            # If no taxons or no data, initialize with empty data
            for taxon in taxons: # Will iterate only if taxons list is not empty
                 divergence_data[taxon] = {
                     'Total score': 0.0,
                     'No of charsets': 0,
                 }
            # Even if taxons is empty, return empty dict
            return divergence_data


        reference_taxon = taxons[0]
        ref_seq = concatenated_sequences.get(reference_taxon, "")
        total_alignment_length = len(ref_seq) # Length of the full concatenated alignment

        # Initialize data structure for each taxon, including N/A for all gene scores initially
        for taxon in taxons:
             divergence_data[taxon] = {
                  'Total score': 0.0,
                  'No of charsets': 0, # Counts genes where this taxon has non-gap data
             }
             # Initialize gene scores as "N/A" using the actual gene names
             for gene in gene_info:
                 divergence_data[taxon][gene['name']] = "N/A"


        # Calculate per-gene differences and count non-gap charsets
        # Only consider genes that have a positive length in the concatenated alignment
        genes_with_length = [gene for gene in gene_info if gene.get('length', 0) > 0]

        for gene in genes_with_length:
             gene_name = gene['name']
             start_0based = gene['start'] - 1 # Convert to 0-based index
             end_0based = gene['end'] # End is inclusive in 1-based, so slice up to end (exclusive)

             if start_0based < 0 or end_0based <= start_0based or end_0based > total_alignment_length:
                  print(f"Warning: Invalid gene range for '{gene_name}' ({gene['start']}-{gene['end']}) within alignment length {total_alignment_length}. Skipping divergence calculation for this gene.", file=sys.stderr)
                  # Mark this gene as N/A for all taxons explicitly if range is bad
                  for taxon in taxons:
                       divergence_data[taxon][gene_name] = "N/A (Bad Range)"
                  continue


             ref_segment = ref_seq[start_0based:end_0based]
             gene_segment_length = len(ref_segment) # Should match gene.get('length', 0)

             if gene_segment_length == 0:
                  # This case should ideally be filtered out by `genes_with_length` list
                  # But as a safeguard, skip if length is zero.
                  print(f"Warning: Gene '{gene_name}' has calculated length 0. Skipping divergence calculation for this gene.", file=sys.stderr)
                  # Mark this gene as N/A for all taxons explicitly
                  for taxon in taxons:
                       divergence_data[taxon][gene_name] = "N/A (Zero Length)"
                  continue

             # Calculate scores for each taxon relative to the reference taxon for this gene segment
             for taxon in taxons:
                  # Reference taxon score is always 0
                  if taxon == reference_taxon:
                       divergence_data[taxon][gene_name] = f"0% #0 (0.00%)"
                       # No of charsets will be calculated *after* the loop, based on which genes have non-gap data

                       continue # Skip calculation for reference taxon


                  taxon_seq = concatenated_sequences.get(taxon, "")
                  # Get the corresponding segment for the current taxon.
                  # Concatenation logic should ensure taxon_seq is long enough and padded.
                  if len(taxon_seq) < end_0based:
                       print(f"Internal Error: Taxon '{taxon}' sequence unexpectedly shorter than end position for gene '{gene_name}'. Cannot calculate divergence.", file=sys.stderr)
                       divergence_data[taxon][gene_name] = "Error (Seq Short)"
                       continue

                  taxon_segment = taxon_seq[start_0based:end_0based]

                  diff_count = 0
                  # comparable_sites = 0 # Not needed for the requested percentage calculation
                  taxon_has_non_gap_in_gene_segment = False

                  # Compare characters in the segment
                  for j in range(gene_segment_length):
                       ref_char = ref_segment[j]
                       taxon_char = taxon_segment[j]

                       # Check if the taxon has *any* non-gap character in this gene segment
                       # This is used to increment the 'No of charsets' count
                       if taxon_char not in ('-', '?'):
                            taxon_has_non_gap_in_gene_segment = True

                       # Only compare if *both* reference and taxon have valid characters at this position
                       if ref_char not in ('-', '?') and taxon_char not in ('-', '?'):
                            # comparable_sites += 1 # This is used for p-distance, but the request is simpler percentage
                            if ref_char.upper() != taxon_char.upper():
                                 diff_count += 1

                  # Calculate percentage based on *gene segment length* (total positions), matching the image format logic
                  # This is NOT pairwise deletion percentage (diff_count / comparable_sites)
                  percentage = (diff_count / gene_segment_length) * 100 if gene_segment_length > 0 else 0.0

                  # Format the score string (e.g., "10% #5 (10.50%)")
                  # The integer part seems to be the first number, followed by count and precise percentage in parentheses
                  leading_percent_int = int(round(percentage)) # Round to nearest integer for the leading part
                  if leading_percent_int > 99: leading_percent_int = 99 # Cap at 99? Image shows max 99. Let's follow.
                  if leading_percent_int < 0: leading_percent_int = 0

                  # Store the formatted gene score using the actual gene name as the key
                  divergence_data[taxon][gene_name] = f"{leading_percent_int}% #{diff_count} ({percentage:.2f}%)"

                  # Increment 'No of charsets' for this taxon if it had non-gap data in this gene
                  if taxon_has_non_gap_in_gene_segment:
                       divergence_data[taxon]['No of charsets'] += 1


        # Calculate Total Score per taxon (Average of gene percentages where taxon had non-gap data)
        for taxon in taxons:
             if taxon == reference_taxon:
                  divergence_data[taxon]['Total score'] = 0.0
                  continue

             total_percentage_sum = 0.0
             genes_counted_for_total_score = 0

             # Iterate through genes that actually had length > 0
             for gene in genes_with_length:
                  gene_name = gene['name']
                  gene_score_str = divergence_data[taxon].get(gene_name)

                  # Check if the score was successfully calculated (not "N/A" or "Error")
                  if gene_score_str and "N/A" not in gene_score_str and "Error" not in gene_score_str:
                       # Extract the precise percentage from the parentheses
                       match = re.search(r'\((\d+\.?\d*)%\)', gene_score_str)
                       if match:
                            try:
                                 percentage_value = float(match.group(1))
                                 total_percentage_sum += percentage_value
                                 genes_counted_for_total_score += 1
                            except ValueError:
                                 print(f"Warning: Could not parse percentage from gene score string for taxon '{taxon}', gene '{gene_name}': {gene_score_str}", file=sys.stderr)
                       else:
                            # Handle cases where the format might be unexpected but not explicitly N/A/Error
                            print(f"Warning: Unexpected gene score format for taxon '{taxon}', gene '{gene_name}': '{gene_score_str}'. Skipping for total score calculation.", file=sys.stderr)


             average_percentage = (total_percentage_sum / genes_counted_for_total_score) if genes_counted_for_total_score > 0 else 0.0
             divergence_data[taxon]['Total score'] = round(average_percentage, 2) # Round final total score


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

    # Raw contents and their assigned names (must match index)
    raw_contents = [
        dummy_fasta_gene1,
        dummy_fasta_gene2,
        dummy_fasta_gene3,
    ]

    assigned_gene_names = [
        "Gene_One",
        "ImportantGene",
        "MyGene_v2",
    ]

    try:
        concatenator = SequenceConcatenator(raw_contents, assigned_gene_names)

        concatenated = concatenator.get_concatenated_sequences()
        stats = concatenator.get_statistics()
        partition = concatenator.get_partition()

        print("\nConcatenated Sequences:")
        if concatenated:
             for taxon, seq in concatenated.items():
                 print(f">{taxon}\n{seq}")
        else:
             print("No concatenated sequences produced.")

        print("\nStatistics:")
        import json
        # Custom JSONEncoder to handle non-serializable types like sets if they appear unexpectedly
        class SetEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, set):
                    return list(obj)
                # Attempt to serialize other potentially complex types
                try:
                    return json.JSONEncoder.default(self, obj)
                except TypeError:
                    return str(obj) # Fallback to string representation

        print(json.dumps(stats, indent=2, cls=SetEncoder))

        print("\nPartition:")
        if partition:
             for gene_name, pos_range, gene_type in partition:
                 print(f"Gene: {gene_name}, Range: {pos_range}, Type: {gene_type}")
        else:
             print("No partition data produced.")

        # Verify divergence data uses correct gene names
        if 'Divergence Data' in stats and stats['Divergence Data']:
             print("\nDivergence Data Check:")
             sample_taxon, data = next(iter(stats['Divergence Data'].items()))
             print(f"Sample Taxon: {sample_taxon}")
             print("Keys in sample data (should include gene names):")
             print(list(data.keys()))
             # Check if the assigned gene names are present as keys (excluding standard keys)
             # Filter assigned_gene_names to include only those that resulted in columns (i.e., genes with length > 0)
             processed_gene_names_with_length = [g['name'] for g in concatenator._SequenceConcatenator__gene_info if g.get('length', 0) > 0] # Access private member for check

             gene_keys_present = sorted([key for key in data.keys() if key not in ['Total score', 'No of charsets']])

             print(f"Gene keys found: {gene_keys_present}")
             print(f"Processed gene names with length > 0: {sorted(processed_gene_names_with_length)}")

             if gene_keys_present == sorted(processed_gene_names_with_length):
                  print("Check Passed: Divergence data keys match processed gene names with length > 0.")
             else:
                  print("Check Failed: Divergence data keys do NOT exactly match processed gene names with length > 0.")
                  print("Missing from divergence keys:", set(processed_gene_names_with_length) - set(gene_keys_present))
                  print("Unexpected in divergence keys:", set(gene_keys_present) - set(processed_gene_names_with_length))

        else:
             print("\nNo Divergence Data produced.")


    except ValueError as e:
        print(f"Configuration Error: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during backend processing:\n{e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)