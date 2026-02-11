#!/usr/bin/env python3
"""Parse a multi-sample VCF file and produce an XML file with DNA sequences."""

import argparse
import re
import sys

from pysam import VariantFile
import numpy as np

datatype = {"binary": "binary", "nd16": "nucleotideDiploid16"}


nd16phased = {
    "AA": "0",
    "AC": "1",
    "AG": "2",
    "AT": "3",
    "CA": "4",
    "CC": "5",
    "CG": "6",
    "CT": "7",
    "GA": "8",
    "GC": "9",
    "GG": "a",
    "GT": "b",
    "TA": "c",
    "TC": "d",
    "TG": "e",
    "TT": "f",
    # Phased with one unknown
    ".A": "g",
    ".C": "h",
    ".G": "i",
    ".T": "j",
    "A.": "n",
    "C.": "o",
    "G.": "p",
    "T.": "q",
    "..": "-",
}


nd16unphased = {
    "AA": "0",
    "CC": "5",
    "GG": "a",
    "TT": "f",
    "AC": "M",
    "CA": "M",
    "AG": "R",
    "GA": "R",
    "AT": "W",
    "TA": "W",
    "CG": "S",
    "GC": "S",
    "CT": "Y",
    "TC": "Y",
    "GT": "K",
    "TG": "K",
    "..": "-",
}


binary = {
    "00": "0",
    "11": "1",
    "01": "1",
    "10": "1",
    ".1": "1",
    "1.": "1",
    "0.": "?",
    ".0": "?",
    "..": "?",
}

ASCERTAINMENT_ADDITION_PHASED = "0123567abf"

ASCERTAINMENT_ADDITION_UNPHASED = "05afMRWSYK"


def parse_args():
    parser = argparse.ArgumentParser(
        "Transform a multi-sample VCF file into an XML file"
    )
    parser.add_argument("vcf", type=str, help="A multi-sample VCF file")
    parser.add_argument("xml_template", type=str, help="An existing XML template file")
    parser.add_argument("output_xml", type=str, help="Path to the output XML file")
    parser.add_argument(
        "--encoding",
        type=str,
        choices=["nd16", "binary"],
        default="nd16",
        help="""Type of encoding for translating diploid variants into a single character.""",
    )
    parser.add_argument(
        "--filtering_density",
        type=float,
        default=0.0,
        help="The required alignment density after filtering",
    )
    return parser.parse_args()


def vcf2xml(vcf, xml_template, output_xml, encoding="nd16", filtering_density=0.0):
    # Parse the VCF file
    names, sequences = parse_vcf(vcf, encoding)
    # Exclude the last sequence --> because its the outgroup

    # names = names[:-1]
    # sequences = sequences[:-1]

    # Read the XML template
    with open(xml_template, "r", encoding="utf-8") as template_file:
        xml_content = template_file.read()

    ascertained = 'ascertained="true"' in xml_content

    seq_length = len(sequences[0])
    print(
        "The initial snv position number was ",
        seq_length,
        "with ",
        count_missing(sequences),
        " missing data.",
    )

    if not 0 <= filtering_density <= 1:
        sys.exit("Error: --filtering_density must be between 0 and 1 (inclusive).")

    if filtering_density != 0:
        # Convert list of strings to a 2D numpy array (shape: n_samples x seq_length)
        arr = np.array([list(seq) for seq in sequences])

        # Compute the proportion of '-' per column
        gap_fraction = np.mean(arr == "-", axis=0)

        # Keep columns where proportion of '-' ≤ density
        keep_mask = gap_fraction <= 1 - filtering_density

        # Rebuild filtered sequences
        sequences = ["".join(row[keep_mask]) for row in arr]

    # Transpose sequences to iterate by columns
    columns = zip(*sequences)

    # Filter columns
    filtered_columns = [
        col
        for col in columns
        if len(set(col) - {"-"}) > 1  # more than one non-gap character
    ]
    # Check if any columns remain after filtering
    if len(filtered_columns) == 0:
        sys.exit("Error: No columns left after filtering. Adjust --filtering_density or check your data.")

    # Rebuild sequences (transpose back)
    sequences = ["".join(chars) for chars in zip(*filtered_columns)]

    seq_length = len(sequences[0])
    print(
        "After filtering ",
        seq_length,
        " snv positions remained with ",
        count_missing(sequences),
        " missing data.",
    )

    if ascertained:
        sequences = [s + ASCERTAINMENT_ADDITION_UNPHASED for s in sequences]

    # Generate new sequence entries
    sequence_entries = []
    for name, seq in zip(names, sequences):
        sequence_entries.append(
            f'\t\t<sequence id="seq_{name}" spec="Sequence" \
                taxon="{name}" totalcount="16" value="{seq}"/>'
        )
    sequences_xml = "\n".join(sequence_entries)

    # Replace the <data id="alignment"> section
    pattern = r'(<data\s+id="alignment"[^>]*>)(.*?)(</data>)'
    if re.search(pattern, xml_content, flags=re.DOTALL):
        # Replace content within <data id="alignment">
        xml_content = re.sub(
            pattern, f"\\1\n{sequences_xml}\n\t\\3", xml_content, flags=re.DOTALL
        )
    else:
        raise ValueError(
            'The XML template must contain a <data id="alignment"> section.'
        )

    if ascertained:
        xml_content = xml_content.replace("insert_from", str(seq_length))
        xml_content = xml_content.replace("insert_to", str(seq_length + 10))

    # Debugging: Print the modified XML content
    # print("Modified XML Content:\n", xml_content)

    # Write the modified XML to the output file
    with open(output_xml, "w", encoding="utf-8") as output_file:
        output_file.write(xml_content)


def parse_vcf(file, encoding):
    with VariantFile(file) as vcf:
        names = list(vcf.header.samples)
        sequences = get_sequences(vcf.fetch(), encoding)
    return names, sequences


def get_sequences(variants, encoding="nd16"):
    sequences = [variant2bases(variant, encoding) for variant in variants]
    sequences = transpose_list(sequences)
    sequences = map("".join, sequences)
    return list(sequences)


def transpose_list(lst):
    return list(map(list, zip(*lst)))


def variant2bases(variant, encoding="nd16"):
    return [sample2base(sample, encoding) for sample in variant.samples.values()]


def sample2base(sample, encoding="nd16"):
    if encoding == "binary":
        alleles = [
            1 if allele is not None and allele > 1 else allele
            for allele in sample.allele_indices
        ]
    else:
        alleles = sample.alleles
    genome = ["." if allele is None else str(allele) for allele in alleles]
    genome = "".join(genome)
    return translate_genome(genome, encoding, False)  # , sample.phased)


def translate_genome(genome, encoding="nd16", phased=True):
    if len(genome) != 2:
        raise ValueError(f"Genome must be diploid. The ploidy is: {len(genome)}")
    if encoding == "binary":
        return binary[genome]
    if encoding == "nd16" and phased:
        return nd16phased[genome]
    if encoding == "nd16" and not phased:
        return nd16unphased[genome]
    else:
        raise ValueError(f'Encoding "{encoding}" is not supported.')


def count_missing(sequences):
    total_gaps = sum(seq.count("-") for seq in sequences)
    return round(total_gaps / (len(sequences) * len(sequences[0])), 3)


if __name__ == "__main__":
    args = parse_args()
    if not 0 <= args.filtering_density <= 1:
        sys.exit("Error: --filtering_density must be between 0 and 1 (inclusive).")
    vcf2xml(
        args.vcf,
        args.xml_template,
        args.output_xml,
        args.encoding,
        args.filtering_density,
    )
