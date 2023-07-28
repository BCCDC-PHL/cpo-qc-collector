#!/usr/bin/env python

import argparse
import datetime
import json
import logging
import os
import time

import cpo_qc_collector.config
import cpo_qc_collector.core as core

DEFAULT_SCAN_INTERVAL_SECONDS = 3600.0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config')
    parser.add_argument('--log-level')
    args = parser.parse_args()

    config = {}
    scan_interval = DEFAULT_SCAN_INTERVAL_SECONDS

    try:
        log_level = getattr(logging, args.log_level.upper())
    except AttributeError as e:
        log_level = logging.INFO

    logging.basicConfig(
        format='{"timestamp": "%(asctime)s.%(msecs)03d", "level": "%(levelname)s", "module", "%(module)s", "function_name": "%(funcName)s", "line_num", %(lineno)d, "message": %(message)s}',
        datefmt='%Y-%m-%dT%H:%M:%S',
        encoding='utf-8',
        level=log_level,
    )
    logging.debug(json.dumps({"event_type": "debug_logging_enabled"}))

    quit_when_safe = False

    while(True):
        try:
            if args.config:
                try:
                    config = cpo_qc_collector.config.load_config(args.config)
                    # Un-comment below to debug config parsing
                    # config['excluded_runs'] = list(config['excluded_runs'])
                    # print(json.dumps(config, indent=2))
                    # exit()
                    logging.info(json.dumps({"event_type": "config_loaded", "config_file": os.path.abspath(args.config)}))
                except json.decoder.JSONDecodeError as e:
                    # If we fail to load the config file, we continue on with the
                    # last valid config that was loaded.
                    logging.error(json.dumps({"event_type": "load_config_failed", "config_file": os.path.abspath(args.config)}))

            core.create_output_dirs(config)

            scan_start_timestamp = datetime.datetime.now()

            runs = core.find_runs(config)
            runs_output_file = os.path.join(config['output_dir'], 'illumina_runs.json')
            with open(runs_output_file, 'w') as f:
                json.dump(runs, f, indent=2)
            logging.info(json.dumps({"event_type": "write_runs_file_complete", "runs_file": runs_output_file}))

            for run in core.scan(config):
                if run is not None:
                    try:
                        config = cpo_qc_collector.config.load_config(args.config)
                        logging.info(json.dumps({"event_type": "config_loaded", "config_file": os.path.abspath(args.config)}))
                    except json.decoder.JSONDecodeError as e:
                        logging.error(json.dumps({"event_type": "load_config_failed", "config_file": os.path.abspath(args.config)}))
                    core.collect_outputs(config, run)
                if quit_when_safe:
                    exit(0)
            scan_complete_timestamp = datetime.datetime.now()
            scan_duration_delta = scan_complete_timestamp - scan_start_timestamp
            scan_duration_seconds = scan_duration_delta.total_seconds()
            logging.info(json.dumps({"event_type": "scan_complete", "scan_duration_seconds": scan_duration_seconds}))

            if quit_when_safe:
                exit(0)

            if "scan_interval_seconds" in config:
                try:
                    config['scan_interval_seconds'] = float(str(config['scan_interval_seconds']))
                except ValueError as e:
                    config['scan_interval_seconds'] = DEFAULT_SCAN_INTERVAL_SECONDS
            else:
                    config['scan_interval_seconds'] = DEFAULT_SCAN_INTERVAL_SECONDS
            time.sleep(config['scan_interval_seconds'])
        except KeyboardInterrupt as e:
            logging.info(json.dumps({"event_type": "quit_when_safe_enabled"}))
            quit_when_safe = True

if __name__ == '__main__':
    main()
