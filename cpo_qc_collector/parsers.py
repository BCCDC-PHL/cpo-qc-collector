import csv

def parse_fastp(fastp_path):
    """
    """
    fastp = []
    int_fields = [
        'total_reads_before_filtering',
        'total_reads_after_filtering',
        'total_bases_before_filtering',
        'total_bases_after_filtering',
        'read1_mean_length_before_filtering',
        'read1_mean_length_after_filtering',
        'read2_mean_length_before_filtering',
        'read2_mean_length_after_filtering',
        'q20_bases_before_filtering',
        'q20_bases_after_filtering',
        'q30_bases_before_filtering',
        'q30_bases_after_filtering',
        'adapter_trimmed_reads',
        'adapter_trimmed_bases',
    ]

    float_fields = [
        'q20_rate_before_filtering',
        'q20_rate_after_filtering',
        'q30_rate_before_filtering',
        'q30_rate_after_filtering',
        'gc_content_before_filtering',
        'gc_content_after_filtering',
    ]
    
    with open(fastp_path, 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            for field in int_fields:
                try:
                    row[field] = int(row[field])
                except ValueError as e:
                    row[field] = None
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    row[field] = None

            fastp.append(row)

    return fastp


def parse_quast(quast_path):
    """
    """
    quast = []

    int_fields = [
        'total_length',
        'num_contigs',
        'largest_contig',
        'assembly_N50',
        'assembly_N75',
        'assembly_L50',
        'assembly_L75',
        'num_contigs_gt_0_bp',
        'num_contigs_gt_1000_bp',
        'num_contigs_gt_5000_bp',
        'num_contigs_gt_10000_bp',
        'num_contigs_gt_25000_bp',
        'num_contigs_gt_50000_bp',
        'total_length_gt_0_bp',
        'total_length_gt_1000_bp',
        'total_length_gt_5000_bp',
        'total_length_gt_10000_bp',
        'total_length_gt_25000_bp',
        'total_length_gt_50000_bp',

    ]
    float_fields = [
        'num_N_per_100_kb',
    ]
    with open(quast_path, 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            for field in int_fields:
                try:
                    row[field] = int(row[field])
                except ValueError as e:
                    row[field] = None
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    row[field] = None

            quast.append(row)

    return quast


def parse_mlst_sequence_type(mlst_sequence_type_path):
    """
    """
    mlst = []
    float_fields = [
        'score',
    ]
    with open(mlst_sequence_type_path, 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    row[field] = None
            mlst.append(row)

    return mlst


def parse_resistance_gene_report(resistance_gene_report_path):
    """
    """
    resistance_genes = []
    int_fields = [
        'resistance_gene_contig_position_start',
        'resistance_gene_contig_position_end',
        'num_contigs_in_plasmid_reconstruction',
        'plasmid_reconstruction_size',
        'depth_coverage_threshold',
        'num_snps_vs_ref_plasmid',
    ]

    float_fields = [
        'percent_resistance_gene_coverage',
        'percent_resistance_gene_identity',
        'mash_neighbor_distance',
        'percent_ref_plasmid_coverage_above_depth_threshold',
    ]
    with open(resistance_gene_report_path, 'r') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        for row in reader:
            for field in int_fields:
                try:
                    row[field] = int(row[field])
                except ValueError as e:
                    row[field] = None
            for field in float_fields:
                try:
                    row[field] = float(row[field])
                except ValueError as e:
                    row[field] = None

            resistance_genes.append(row)

    return resistance_genes
                
            
