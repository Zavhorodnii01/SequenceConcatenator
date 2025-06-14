# Import the regular expression module for pattern matching in sequences and headers.
import re

# Define the main class responsible for concatenating gene sequences from various formats.
class SequenceConcatenator:
    """
    A class to parse gene sequence files (FASTA, Nexus, GenBank),
    concatenate sequences for common taxons, and generate statistics
    and partition information.
    """

    # Constructor method to initialize the SequenceConcatenator.
    def __init__(self, gene_file_contents: list[list[str]]):
        """
        Initializes the SequenceConcatenator with raw gene file contents.

        Args:
            gene_file_contents: A list where each element is a list of strings,
                                representing the lines of a single gene file.
        """

        self.__raw_gene_contents = gene_file_contents
        self.__parsed_gene_data = []
        self.__gene_info = []
        self.__all_taxons = set()
        self.__concatenated_sequences = {}
        self.__partition_data = []
        self.__statistics = {}

        self._parse_all_genes()
        self._collect_all_taxons()
        self.__concatenated_sequences = self._concatenate_sequences()
        self.__partition_data = self._calculate_partition()
        self.__statistics = self._calculate_statistics()

    def get_concatenated_sequences(self) -> dict[str, str]:
        """
        Returns the dictionary of concatenated sequences.

        Returns:
            A dictionary where keys are taxon names and values are their
            concatenated sequences.
        """
        return self.__concatenated_sequences

    def get_statistics(self) -> dict:
        """
        Returns the calculated statistics for the concatenated alignment.

        Returns:
            A dictionary containing various statistics like number of taxa,
            total length, missing data, gene lengths, etc.
        """
        return self.__statistics

    def get_partition(self) -> list[tuple[str, str, str]]:
        """
        Returns the partition data for the concatenated alignment.

        Returns:
            A list of tuples, where each tuple contains (gene_name, position_range, gene_type).
        """
        return self.__partition_data

    def _parse_all_genes(self) -> None:
        """
        Iterates through all raw gene file contents, detects the file format
        (FASTA, Nexus, GenBank), and calls the appropriate parsing method.
        Stores the parsed data and basic gene info.
        """
        # Iterates through each gene file's content and its index.
        for i, file_lines in enumerate(self.__raw_gene_contents):
            content_str = "".join(file_lines).strip()

            gene_name = f"gene{i + 1}"
            try:
                # Attempts to detect the file format based on common headers/keywords.
                if content_str.lower().startswith('#nexus') or 'begin data' in content_str.lower():
                    parsed_data = self._parse_nexus(file_lines)
                elif content_str.lower().startswith('locus') or content_str.lower().startswith(
                        'version') or 'ncbi' in content_str.lower():
                    parsed_data = self._parse_genbank(file_lines)
                elif content_str.startswith('>'):
                    parsed_data = self._parse_fasta(file_lines)
                else: # Otherwise
                    parsed_data = self._parse_fasta(file_lines)
                if not parsed_data:
                    continue

                # Finds first sequence in parsed_data
                first_seq = ""
                if parsed_data:
                     first_seq = list(parsed_data.values())[0]

                seq_chars_str = "".join(first_seq).replace('-', '').replace('?', '').replace(' ', '').replace('\t','').strip().upper()

                # Defines sets of standard DNA and Protein characters.
                dna_chars = set("ACGTU")
                protein_chars = set("ACDEFGHIKLMNPQRSTVWYBZX")

                # is Used to ensure data input is correct
                other_chars = set(seq_chars_str) - (dna_chars | protein_chars)

                # Determines sequence type based on characters present.
                if len(seq_chars_str) > 0 and all(c in dna_chars for c in seq_chars_str):
                    determined_seq_type = 'DNA'
                elif len(seq_chars_str) > 0 and not other_chars:
                    protein_overlap = protein_chars & set(seq_chars_str)
                    if 'T' in seq_chars_str or 'U' in seq_chars_str or (len(protein_overlap) > 0 and len(protein_overlap) < len(set(seq_chars_str))/2):
                        determined_seq_type = 'DNA'
                    else:
                        determined_seq_type = 'Protein'
                else:
                    determined_seq_type = 'Unknown'

                self.__parsed_gene_data.append(parsed_data)
                self.__gene_info.append({'name': gene_name, 'type': determined_seq_type})

            except Exception as e:
                pass

    def _parse_fasta(self, lines: list[str]) -> dict[str, str]:
        """
        Parses sequence data from a list of lines assuming FASTA format.

        Args:
            lines: A list of strings representing the lines of a FASTA file.

        Returns:
            A dictionary where keys are sequence names (from headers) and
            values are the corresponding concatenated sequences. Returns an
            empty dictionary if no valid FASTA data is found.
        """
        sequences = {}
        current_name = None
        current_seq_lines = []

        lines = [line.strip() for line in lines]

        for line in lines:
            if not line or line.startswith('#'):
                continue
            if line.startswith(">"):
                if current_name is not None and current_seq_lines:
                    sequences[current_name] = "".join(current_seq_lines)
                match = re.match(r'>\s*(.+)', line)
                if match:
                    full_header = match.group(1).strip()

                    current_name = full_header
                    if not current_name:
                        current_name = None
                        current_seq_lines = []
                        continue
                else:
                    current_name = line[1:].strip()
                    if not current_name:
                        current_name = None
                        current_seq_lines = []
                        continue
                current_seq_lines = []
            elif current_name is not None:
                clean_seq_line = re.sub(r'[\s\d]+', '', line)
                current_seq_lines.append(clean_seq_line)

        if current_name is not None and current_seq_lines:
            sequences[current_name] = "".join(current_seq_lines)

        return sequences

    def _parse_nexus(self, lines: list[str]) -> dict[str, str]:
        """
        Parses sequence data from a list of lines assuming Nexus format,
        specifically looking for the 'MATRIX' block.

        Args:
            lines: A list of strings representing the lines of a Nexus file.

        Returns:
            A dictionary where keys are taxon names and values are their
            corresponding sequences from the MATRIX block. Returns an
            empty dictionary if no valid MATRIX block is found.
        """
        sequences = {}
        in_matrix = False
        current_taxon = None
        current_seq_parts = []

        for line in lines:
            line = line.strip()
            if line.startswith('#') or line.startswith('[') or not line:
                continue

            line_lower = line.lower()
            if line_lower.startswith('begin data;'):

                sequences = {}
                in_matrix = False
                continue


            if line_lower.startswith('matrix;'):
                in_matrix = True
                current_taxon = None
                current_seq_parts = []
                continue

            if line_lower.startswith('end;'):
                if in_matrix and current_taxon is not None and current_seq_parts:
                    sequences[current_taxon] = "".join(current_seq_parts)
                    current_taxon = None

                if in_matrix:
                    break
                continue
            if in_matrix:
                parts = line.split()
                if not parts: continue

                first_part = parts[0]

                name_match = re.match(r"^(['\"]).*?\1", first_part)
                if name_match:
                    taxon_name = first_part.strip(name_match.group(1))
                    sequence_part_raw = "".join(parts[1:])
                else:
                    taxon_name = first_part
                    sequence_part_raw = "".join(parts[1:])

                if taxon_name:
                    sequence_part = sequence_part_raw.replace(';', '').strip()

                    if current_taxon is None or current_taxon != taxon_name:
                        if current_taxon is not None and current_seq_parts:
                            sequences[current_taxon] = "".join(current_seq_parts)
                        current_taxon = taxon_name
                        current_seq_parts = [sequence_part]
                    elif current_taxon == taxon_name:
                         current_seq_parts.append(sequence_part)

        if current_taxon is not None and current_seq_parts:
            sequences[current_taxon] = "".join(current_seq_parts)

        return sequences

    def _parse_genbank(self, lines: list[str]) -> dict[str, str]:
        """
        Parses sequence data from a list of lines assuming GenBank format.
        Extracts the organism name and the sequence from the ORIGIN block.

        Args:
            lines: A list of strings representing the lines of a GenBank file.

        Returns:
            A dictionary where the key is the organism name and the value
            is the full sequence. Returns an empty dictionary if no organism
            or origin block with sequence data is found.
        """
        sequences = {}
        current_organism_name = None
        in_origin = False
        sequence_lines = []


        for i, line in enumerate(lines):
            line_stripped = line.rstrip()
            if not line_stripped: continue

            if line_stripped == '//':
                if current_organism_name is not None and sequence_lines:
                    sequence = "".join(sequence_lines).replace(' ', '').replace('\t', '').replace('\n', '')

                    if sequence:
                         sequences[current_organism_name] = sequence


                current_organism_name = None
                in_origin = False
                sequence_lines = []
                continue


            org_match = re.match(r'^\s{2,}ORGANISM\s+(.*?)\s*\.?$', line, re.IGNORECASE)
            if org_match:
                current_organism_name = org_match.group(1).strip()
                in_origin = False
                sequence_lines = []
                continue


            if line_stripped.lower() == 'origin':
                in_origin = True
                sequence_lines = []
                continue


            if in_origin:

                seq_part = re.sub(r'^\s*\d+\s*', '', line).replace(' ', '').replace('\t', '')
                sequence_lines.append(seq_part)
                continue
        if current_organism_name is not None and sequence_lines and (not lines or '//' not in lines[-1][-5:]):
             sequence = "".join(sequence_lines).replace(' ', '').replace('\t', '').replace('\n', '')
             if sequence:
                 if current_organism_name in sequences:

                    sequences[current_organism_name] = sequence
                 else:

                     sequences[current_organism_name] = sequence
        return sequences

    def _collect_all_taxons(self) -> None:
        """
        Collects all unique taxon names from the parsed data of all genes
        and stores them in a sorted list.
        """
        self.__all_taxons = set()
        for gene_data_dict in self.__parsed_gene_data:
            self.__all_taxons.update(gene_data_dict.keys())
        self.__all_taxons = sorted(list(self.__all_taxons))

    def _concatenate_sequences(self) -> dict[str, str]:
        """
        Concatenates the sequences for each taxon across all parsed genes.
        If a taxon is missing from a gene, its sequence for that gene's
        length is filled with hyphens ('-'). Updates gene info with
        start/end positions in the concatenated alignment.

        Returns:
            A dictionary where keys are all unique taxon names and values
            are their concatenated sequences.
        """
        concatenated_sequences = {}

        # Uses dictionary to write taxons names as keys
        for taxon in self.__all_taxons:
            concatenated_sequences[taxon] = ""

        current_position = 1

        # Iterates through the parsed data for each gene.
        for i, gene_data_dict in enumerate(self.__parsed_gene_data):
            gene_length = 0

            parsed_data = gene_data_dict
            if parsed_data:
                # Gets the first taxon name
                first_taxon_in_gene = list(parsed_data.keys())[0]

                # Gets the length of that taxon's sequence.
                gene_length = len(parsed_data[first_taxon_in_gene])

            # Iterates through all unique taxons identified across all genes.
            for taxon in self.__all_taxons:
                # Takes the current sequence from current taxon in current gene
                sequence = gene_data_dict.get(taxon)

                if sequence is not None: # The length check
                    if len(sequence) != gene_length and gene_length > 0:
                        # Inserting missing data
                         if len(sequence) < gene_length:
                             sequence += "-" * (gene_length - len(sequence))
                         # The case where the sequence is longer than expected
                         elif len(sequence) > gene_length:
                             sequence = sequence[:gene_length]
                         concatenated_sequences[taxon] += sequence

                    # if the sequence is present and its full, just inserting the whole sequence
                    else:
                         concatenated_sequences[taxon] += sequence
                # if taxon isn't found
                else:
                    # Missing data -
                    concatenated_sequences[taxon] += "-" * gene_length

            # Updates gene's info about lenght , start and end
            if i < len(self.__gene_info):
                self.__gene_info[i]['length'] = gene_length
                self.__gene_info[i]['start'] = current_position
                self.__gene_info[i]['end'] = current_position + gene_length - 1
                # Updates the current_position for the next gene.
                current_position += gene_length
        return concatenated_sequences

    def _calculate_statistics(self) -> dict:
        """
        Calculates various statistics about the concatenated sequence alignment.

        Returns:
            A dictionary containing statistics such as number of taxa, total length,
            missing data counts, percentage missing data, number of genes,
            and gene lengths within the concatenated alignment.
        """
        statistics = {}

        num_taxa = len(self.__all_taxons)
        statistics["Number of Taxa"] = num_taxa

        total_missing_chars = 0
        missing_data_per_taxon = {}

        if num_taxa > 0:
            arbitrary_taxon = self.__all_taxons[0]
            total_length = len(self.__concatenated_sequences[arbitrary_taxon])
            statistics["Total Length"] = total_length

            for taxon in self.__all_taxons:
                seq = self.__concatenated_sequences[taxon]
                missing_count = seq.count("-") + seq.count("?")
                missing_data_per_taxon[taxon] = missing_count
                total_missing_chars += missing_count

            statistics["Missing Data per Taxon (count)"] = missing_data_per_taxon
            if total_length * num_taxa > 0:
                percentage_missing_overall = (total_missing_chars / (total_length * num_taxa)) * 100
                statistics["Percentage Overall Missing Data (%)"] = round(percentage_missing_overall, 2)
            else:
                statistics["Percentage Overall Missing Data (%)"] = 0.0
        gene_lengths_in_concat = {info['name']: info.get('length', 0) for info in self.__gene_info}
        statistics["Gene Lengths (in concatenated alignment)"] = gene_lengths_in_concat
        statistics["Number of Genes"] = len(self.__gene_info)

        return statistics

    def _calculate_partition(self) -> list[tuple[str, str, str]]:
        """
        Generates partition information for the concatenated alignment
        based on the processed gene information. This is useful for
        phylogenetic software (e.g., RAxML, IQ-TREE).

        Returns:
            A list of tuples, where each tuple represents a partition:
            (gene_name, position_range, gene_type). Positions are 1-based.
        """
        partition = []

        for gene_info in self.__gene_info:
            gene_name = gene_info.get('name', 'UnknownGene')
            start = gene_info.get('start')
            end = gene_info.get('end')
            gene_type = gene_info.get('type', 'Unknown')

            if start is not None and end is not None:
                partition.append((gene_name, f"{start}-{end}", gene_type))
        return partition