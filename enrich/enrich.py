from __future__ import print_function
import argparse
import logging
import json
import config_check
from enrich_error import EnrichError
from experiment import Experiment
from selection import Selection
from seqlib.basic import BasicSeqLib
from seqlib.barcodevariant import BarcodeVariantSeqLib
from seqlib.barcode import BarcodeSeqLib
from seqlib.overlap import OverlapSeqLib

_DRIVER_NAME = "enrich.py"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument("--log", metavar="file", help="path to log file")
    parser.add_argument("--report-filtered-reads", action="store_true", default=False, dest="report_filtered", help="output filtered reads to log file")
    parser.add_argument("--no-plots", help="don't make plots", dest="plots", action="store_false", default=True)
    args = parser.parse_args()
 
    if args.report_filtered:
        log_level = logging.DEBUG
        if args.log is None:
            raise EnrichError("Cannot report filtered reads without a log file", _DRIVER_NAME)
    else:
        log_level = logging.INFO

    if args.log:
        logging.basicConfig(filename=args.log, level=log_level)
    else:
        logging.basicConfig(level=log_level)

    try:
        config = json.load(open(args.config, "U"))
    except IOError:
        raise EnrichError('Failed to open "%s"' % args.config, _DRIVER_NAME)
    except ValueError:
        raise EnrichError("Improperly formatted .json file", _DRIVER_NAME)

    if config_check.is_experiment(config):
        obj = Experiment(config)
    elif config_check.is_selection(config):
        obj = Selection(config)
    elif config_check.is_seqlib(config):
        obj = globals()[config_check.seqlib_type(config)](config)
    else:
        raise EnrichError("Unrecognized .json config", _DRIVER_NAME)

    if obj.output_base is None:
        raise EnrichError("No output directory set", _DRIVER_NAME)

    obj.calculate()
    obj.write_all()
    if args.plots:
        pass
#        obj.make_plots()
    