import collections
import csv
import glob
import json
import logging
import os
import re
import shutil

from pathlib import Path
from typing import Iterator, Optional

import cpo_qc_collector.parsers as parsers


def create_output_dirs(config):
    """
    Create output directories if they don't exist.

    :param config: Application config.
    :type config: dict[str, object]
    :return: None
    :rtype: None
    """
    base_outdir = config['output_dir']
    output_dirs = [
        base_outdir,
        os.path.join(base_outdir, 'library-qc'),
        os.path.join(base_outdir, 'plasmid-qc'),
    ]
    for output_dir in output_dirs:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)    


def find_latest_plasmid_screen_output(analysis_dir: Path, analysis_type: str='short'):
    """
    Find the latest plasmid screen output directory, within the analysis directory.

    :param analysis_dir: Analysis directory.
    :type analysis_dir: str
    :return: Path to latest plasmid screen output directory.
    :rtype: Path
    """
    plasmid_screen_output_dir_glob = "plasmid-screen-v*-output"
    plasmid_screen_output_dirs = glob.glob(os.path.join(analysis_dir, analysis_type, plasmid_screen_output_dir_glob))
    latest_plasmid_screen_output_dir = None
    if len(plasmid_screen_output_dirs) > 0:
        latest_plasmid_screen_output_dir = Path(os.path.abspath(plasmid_screen_output_dirs[-1]))

    return latest_plasmid_screen_output_dir


def find_latest_assembly_output(analysis_dir: Path, analysis_type: str='short'):
    """
    Find the latest routine assembly output directory, within the analysis directory.

    :param analysis_dir: Analysis directory.
    :type analysis_dir: str
    :return: Path to latest routine assembly output directory.
    :rtype: str
    """
    if analysis_type == 'short':
        assembly_output_dir_glob = "routine-assembly-v*-output"
    elif analysis_type == 'hybrid':
        assembly_output_dir_glob = "dragonflye-nf-v*-output"
    assembly_output_dirs = glob.glob(os.path.join(analysis_dir, analysis_type, assembly_output_dir_glob))
    latest_assembly_output_dir = None
    if len(assembly_output_dirs) > 0:
        latest_assembly_output_dir = os.path.abspath(assembly_output_dirs[-1])

    return latest_assembly_output_dir


def find_latest_mlst_output(analysis_dir, analysis_type='short'):
    """
    Find the latest mlst output directory, within the analysis directory.

    :param analysis_dir: Analysis directory.
    :type analysis_dir: str
    :param analysis_type: Analysis type. One of: ['short', 'hybrid', 'long']
    :type analysis_type: str
    :return: Path to latest mlst output directory.
    :rtype: str
    """
    mlst_output_dir_glob = "mlst-nf-v*-output"
    mlst_output_dirs = glob.glob(os.path.join(analysis_dir, analysis_type, mlst_output_dir_glob))
    latest_mlst_output_dir = None
    if len(mlst_output_dirs) > 0:
        latest_mlst_output_dir = os.path.abspath(mlst_output_dirs[-1])

    return latest_mlst_output_dir


def find_analysis_dirs(config):
    """
    Find all analysis directories with completed plasmid-screen analyses.

    :param config: Application config.
    :type config: dict[str, object]
    :return: List of analysis directories.
    :rtype: list[str]
    """
    miseq_run_id_regex = "\d{6}_M\d{5}_\d+_\d{9}-[A-Z0-9]{5}"
    nextseq_run_id_regex = "\d{6}_VH\d{5}_\d+_[A-Z0-9]{9}"
    nanopore_run_id_regex = "\d{8}_\d{4}_X\d+_[A-Z0-9]{8}_[a-z0-9]{8}"

    analysis_by_run_dir = config['analysis_by_run_dir']
    subdirs = os.scandir(analysis_by_run_dir)

    for subdir in subdirs:
        run_id = subdir.name
        matches_miseq_regex = re.match(miseq_run_id_regex, run_id)
        matches_nextseq_regex = re.match(nextseq_run_id_regex, run_id)
        matches_nanopore_regex = re.match(nanopore_run_id_regex, run_id)
        sequencer_type = None
        if matches_miseq_regex:
            sequencer_type = 'miseq'
        elif matches_nextseq_regex:
            sequencer_type = 'nextseq'
        elif matches_nanopore_regex:
            sequencer_type = 'nanopore'
        not_excluded = run_id not in config['excluded_runs']
        ready_to_collect = False

        analysis_types = ['short', 'hybrid', 'long']
        for analysis_type in analysis_types:
            latest_plasmid_screen_output = find_latest_plasmid_screen_output(subdir, analysis_type)
            ready_to_collect = False
            if latest_plasmid_screen_output is not None and os.path.exists(latest_plasmid_screen_output):
                plasmid_screen_complete = os.path.exists(os.path.join(latest_plasmid_screen_output, 'analysis_complete.json'))
                ready_to_collect = plasmid_screen_complete

            conditions_checked = {
                "is_directory": subdir.is_dir(),
                "matches_recognized_run_id_format": ((matches_miseq_regex is not None) or (matches_nextseq_regex is not None) or (matches_nanopore_regex is not None)),
                "not_excluded": not_excluded,
                "ready_to_collect": ready_to_collect,
            }
            conditions_met = list(conditions_checked.values())

            analysis_directory_path = os.path.abspath(subdir.path)
            analysis_dir = {
                "path": analysis_directory_path,
                "sequencer_type": sequencer_type,
                "analysis_type": analysis_type,
            }
            if all(conditions_met):
                logging.info(json.dumps({
                    "event_type": "analysis_directory_found",
                    "sequencing_run_id": run_id,
                    "analysis_type": analysis_type,
                    "analysis_directory_path": analysis_directory_path
                }))

                yield analysis_dir
            else:
                logging.debug(json.dumps({
                    "event_type": "directory_skipped",
                    "analysis_directory_path": os.path.abspath(subdir.path),
                    "conditions_checked": conditions_checked
                }))
                yield None


def find_runs(config):
    """
    Finda all runs that have routine sequence QC data.

    :param config: Application config.
    :type config: dict[str, object]
    :return: List of runs. Keys: ['run_id', 'sequencer_type']
    :rtype: list[dict[str, str]]
    """
    logging.info(json.dumps({"event_type": "find_runs_start"}))
    runs = []
    all_analysis_dirs = sorted(list(os.listdir(config['analysis_by_run_dir'])))
    illumina_run_ids = filter(lambda x: re.match('\d{6}_[VM]', x) != None, all_analysis_dirs)
    nanopore_run_ids = filter(lambda x: re.match('\d{8}_\d{4}_X\d{1}_', x) != None, all_analysis_dirs)
    all_run_ids = sorted(list(illumina_run_ids) + list(nanopore_run_ids))
    for run_id in all_run_ids:
        if run_id in config['excluded_runs']:
            continue

        sequencer_type = None
        if re.match('\d{6}_M\d{5}_', run_id):
            sequencer_type = 'miseq'
        elif re.match('\d{6}_VH\d{5}_', run_id):
            sequencer_type = 'nextseq'
        elif re.match('\d{8}_\d{4}_X\d{1}_', run_id):
            sequencer_type = 'nanopore'

        analysis_dir = os.path.join(config['analysis_by_run_dir'], run_id)

        analysis_types = ['short', 'hybrid', 'long']
        for analysis_type in analysis_types:
            latest_plasmid_screen_output_dir = find_latest_plasmid_screen_output(analysis_dir, analysis_type)

            if latest_plasmid_screen_output_dir is not None and os.path.exists(os.path.join(latest_plasmid_screen_output_dir, 'analysis_complete.json')):
                run = {
                    'run_id': run_id,
                    'sequencer_type': sequencer_type,
                    'analysis_type': analysis_type,
                }
                runs.append(run)

    logging.info(json.dumps({
        "event_type": "find_runs_complete"
    }))

    return runs


def scan(config: dict[str, object]) -> Iterator[Optional[dict[str, str]]]:
    """
    Scanning involves looking for all existing runs and...

    :param config: Application config.
    :type config: dict[str, object]
    :return: A run directory to analyze, or None
    :rtype: Iterator[Optional[dict[str, object]]]
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    for analysis_dir in find_analysis_dirs(config):    
        yield analysis_dir


def infer_species(config: dict[str, object], species_abundance, project_id):
    """
    :return: Name of inferred species.
    :rtype: str
    """
    species_abundance_keys = [
        'abundance_5_',
        'abundance_4_',
        'abundance_3_',
        'abundance_2_',
        'abundance_1_',
    ]
    inferred_species = None
    if 'projects' in config and project_id in config['projects']:
        project = config['projects'][project_id]
        if 'fixed_genome_size' in project and project['fixed_genome_size']:
            inferred_species = project.get('project_species_name', None)
        else:
            greatest_fraction_total_reads = 0.0
            for k in species_abundance_keys:
                if species_abundance[k + 'name'] != 'Homo sapiens' and species_abundance[k + 'fraction_total_reads'] > greatest_fraction_total_reads:
                    inferred_species = species_abundance[k + 'name']
                    greatest_fraction_total_reads = species_abundance[k + 'fraction_total_reads']
    else:
        greatest_fraction_total_reads = 0.0
        for k in species_abundance_keys:
            if species_abundance[k + 'name'] != 'Homo sapiens' and species_abundance[k + 'fraction_total_reads'] > greatest_fraction_total_reads:
                inferred_species = species_abundance[k + 'name']
                greatest_fraction_total_reads = species_abundance[k + 'fraction_total_reads']

    return inferred_species


def get_percent_reads_by_species_name(species_abundance, species_name):
    """
    """
    percent_reads = None
    species_abundance_keys = [
        'abundance_1_',
        'abundance_2_',
        'abundance_3_',
        'abundance_4_',
        'abundance_5_',
    ]

    for k in species_abundance_keys:
        try:
            if species_abundance[k + 'name'] == species_name:
                percent_reads = 100 * species_abundance[k + 'fraction_total_reads']
        except KeyError as e:
            pass

    return percent_reads        
    
    
def collect_outputs(config: dict[str, object], analysis_dir: Optional[dict[str, str]]):
    """
    Collect QC outputs for a specific analysis dir.

    :param config: Application config.
    :type config: dict[str, object]
    :param analysis_dir: Analysis dir. Keys: ['path', 'sequencer_type', 'analysis_type']
    :type analysis_dir: dict[str, str]
    :return: 
    :rtype: 
    """
    run_id = os.path.basename(analysis_dir['path'])
    analysis_type = analysis_dir['analysis_type']
    if analysis_type not in ['short', 'hybrid', 'long']:
        return None
    logging.info(json.dumps({"event_type": "collect_outputs_start", "sequencing_run_id": run_id, "analysis_dir_path": os.path.join(analysis_dir['path'], analysis_type)}))

    # library-qc
    library_qc_by_library_id = {}
    library_qc_dst_file = os.path.join(config['output_dir'], "library-qc", run_id + "_" + analysis_type + "_library_qc.json")
    if not os.path.exists(library_qc_dst_file):
        latest_assembly_output_path = find_latest_assembly_output(analysis_dir['path'], analysis_type)
        if latest_assembly_output_path is not None:
            fastp_glob = os.path.join(latest_assembly_output_path, '*', '*_fastp.csv')
        else:
            logging.error(json.dumps({"event_type": "latest_assembly_output_not_found", "sequencing_run_id": run_id, "analysis_dir_path": os.path.join(analysis_dir['path'], analysis_type)}))
            # if we can't find an assembly dir, bail out early.
            return None
        fastp_paths = glob.glob(fastp_glob)
        for fastp_path in fastp_paths:
            if os.path.exists(fastp_path):
                fastp = parsers.parse_fastp(fastp_path)
                for fastp_record in fastp:
                    library_id = fastp_record['sample_id']
                    if library_id not in library_qc_by_library_id:
                        library_qc_by_library_id[library_id] = {
                            'library_id': library_id,
                            'assembly_type': analysis_type,
                        }
                    library_qc_by_library_id[library_id]['num_bases_short'] = fastp_record['total_bases_before_filtering']

        
        nanoq_glob = os.path.join(latest_assembly_output_path, '*', '*_nanoq.csv')
        nanoq_paths = glob.glob(nanoq_glob)
        for nanoq_path in nanoq_paths:
            if os.path.exists(nanoq_path):
                nanoq = parsers.parse_nanoq(nanoq_path)
                for nanoq_record in nanoq:
                    library_id = nanoq_record['sample_id']
                    if library_id not in library_qc_by_library_id:
                        library_qc_by_library_id[library_id] = {
                            'library_id': library_id,
                            'assembly_type': analysis_type,
                        }
                    library_qc_by_library_id[library_id]['num_bases_long'] = nanoq_record['total_bases_before_filtering']

        quast_glob = os.path.join(latest_assembly_output_path, '*', '*_quast.csv')
        quast_paths = glob.glob(quast_glob)
        for quast_path in quast_paths:
            if os.path.exists(quast_path):
                quast = parsers.parse_quast(quast_path)
                for quast_record in quast:
                    assembly_id_split = quast_record['assembly_id'].split('_')
                    if len(assembly_id_split) > 0:
                        library_id = assembly_id_split[0]
                        if library_id in library_qc_by_library_id:
                            library_qc_by_library_id[library_id]['assembly_length'] = quast_record['total_length']
                            library_qc_by_library_id[library_id]['assembly_num_contigs'] = quast_record['num_contigs']
                            library_qc_by_library_id[library_id]['assembly_N50'] = quast_record['assembly_N50']

        latest_mlst_output_path = find_latest_mlst_output(analysis_dir['path'])
        if latest_mlst_output_path is not None:
            mlst_sequence_type_glob = os.path.join(latest_mlst_output_path, '*', '*_sequence_type.csv')
        else:
            logging.error(json.dumps({"event_type": "latest_mlst_output_not_found", "sequencing_run_id": run_id, "analysis_dir_path": os.path.join(analysis_dir['path'], analysis_type)}))
            mlst_sequence_type_glob = ''
        mlst_sequence_type_paths = glob.glob(mlst_sequence_type_glob)
        for mlst_sequence_type_path in mlst_sequence_type_paths:
            if os.path.exists(mlst_sequence_type_path):
                mlst_sequence_type = parsers.parse_mlst_sequence_type(mlst_sequence_type_path)
                for mlst_record in mlst_sequence_type:
                    library_id = mlst_record['sample_id']
                    if library_id in library_qc_by_library_id:
                        library_qc_by_library_id[library_id]['mlst_scheme'] = mlst_record['scheme']
                        library_qc_by_library_id[library_id]['mlst_sequence_type'] = mlst_record['sequence_type']
                        library_qc_by_library_id[library_id]['mlst_score'] = mlst_record['score']
        
        with open(library_qc_dst_file, 'w') as f:
            json.dump(list(library_qc_by_library_id.values()), f, indent=2)
            f.write('\n')

        logging.info(json.dumps({
            "event_type": "write_library_qc_complete",
            "run_id": run_id,
            "dst_file": library_qc_dst_file
        }))

    # plasmid-qc
    plasmid_qc_by_library_id = {}
    plasmid_qc_dst_file = os.path.join(config['output_dir'], "plasmid-qc", run_id + "_" + analysis_type + "_plasmid_qc.json")
    if not os.path.exists(plasmid_qc_dst_file):
        latest_plasmid_screen_output_path = find_latest_plasmid_screen_output(analysis_dir['path'], analysis_type)
        resistance_gene_glob = os.path.join(latest_plasmid_screen_output_path, '*', '*_resistance_gene_report.tsv')
        resistance_gene_paths = glob.glob(resistance_gene_glob)
        for resistance_gene_path in resistance_gene_paths:
            if os.path.exists(resistance_gene_path):
                resistance_gene_report = parsers.parse_resistance_gene_report(resistance_gene_path)
                for resistance_gene_record in resistance_gene_report:
                    library_id = resistance_gene_record['sample_id']
                    if library_id not in plasmid_qc_by_library_id:
                        plasmid_qc_by_library_id[library_id] = set()
                    resistance_gene = {
                        'library_id': library_id,
                        'assembly_type': analysis_type,
                    }
                    resistance_gene['assembly_type'] = analysis_type
                    resistance_gene['resistance_gene_name'] = resistance_gene_record['resistance_gene_id']
                    resistance_gene['resistance_gene_identity'] = resistance_gene_record['percent_resistance_gene_identity']
                    resistance_gene['mob_suite_primary_cluster_id'] = resistance_gene_record['mob_suite_primary_cluster_id']
                    resistance_gene['mob_suite_secondary_cluster_id'] = resistance_gene_record['mob_suite_secondary_cluster_id']
                    resistance_gene['plasmid_reconstruction_size'] = resistance_gene_record['plasmid_reconstruction_size']
                    resistance_gene['num_contigs_in_plasmid_reconstruction'] = resistance_gene_record['num_contigs_in_plasmid_reconstruction']
                    resistance_gene['closest_db_plasmid'] = resistance_gene_record['mash_nearest_neighbor']
                    
                    plasmid_qc_by_library_id[library_id].add(json.dumps(resistance_gene, sort_keys=True))

        # There are unfortunately sometimes duplicate entries in resistance gene reports
        # We eliminate duplicates by using a set and then converting back to lists.
        for library_id in plasmid_qc_by_library_id.keys():
            resistance_gene_set_of_strings = plasmid_qc_by_library_id[library_id]
            plasmid_qc_by_library_id[library_id] = []
            for resistance_gene_string in resistance_gene_set_of_strings:
                plasmid_qc_by_library_id[library_id].append(json.loads(resistance_gene_string))

        with open(plasmid_qc_dst_file, 'w') as f:
            json.dump(list(plasmid_qc_by_library_id.values()), f, indent=2)
            f.write('\n')

        logging.info(json.dumps({
            "event_type": "write_plasmid_qc_complete",
            "run_id": run_id,
            "dst_file": plasmid_qc_dst_file
        }))
        

    logging.info(json.dumps({"event_type": "collect_outputs_complete", "sequencing_run_id": run_id, "analysis_dir_path": analysis_dir['path']}))
