from __future__ import print_function
import time
import logging
from enrich_error import EnrichError
from datacontainer import DataContainer
import os.path


class SeqLib(DataContainer):
    """
    Abstract class for handling count data from a single sequencing library.
    Creating a :py:class:`seqlib.seqlib.SeqLib` requires a valid *config* 
    object, usually from a ``.json`` configuration file.

    .. note:: Example configuration files can be found in the documentation \
    for derived classes.
    """
    def __init__(self, config):
        DataContainer.__init__(self, config)

        try:
            self.timepoint = int(config['timepoint'])
        except KeyError as key:
            raise EnrichError("Missing required config value '{key}'".format(key=key), 
                              self.name)
        except ValueError as value:
            raise EnrichError("Invalid parameter value {value}".format(value=value), self.name)

        if 'align variants' in config:
            if config['align variants']:
                self.aligner = Aligner()
            else:
                self.aligner = None
        else:
            self.aligner = None

        if 'report filtered reads' in config:
            self.report_filtered_reads = config['report filtered reads']
        else:
            self.report_filtered_reads = self.verbose

        # initialize data
        self.df_dict = dict()        # pandas dataframes
        self.df_file = dict()   # paths to saved counts
        self.filters = None         # dictionary
        self.filter_stats = None    # dictionary


    def calculate(self):
        """
        Pure virtual method that defines how the data are counted.
        """
        raise NotImplementedError("must be implemented by subclass")


    def report_filtered_read(self, fq, filter_flags):
        """
        Write the :py:class:`~fqread.FQRead` object *fq* to the ``DEBUG``
        logging . The dictionary *filter_flags* contains ``True`` 
        values for each filtering option that applies to *fq*. Keys in 
        *filter_flags* are converted to messages using the 
        ``DataContainer._filter_messages`` dictionary.
        """
        logging.debug("Filtered read ({messages}) [{name}]\n{read!s}".format(
                      messages=', '.join(DataContainer._filter_messages[x] 
                                for x in filter_flags if filter_flags[x]), 
                      name=self.name, read=fq))


    def write_all(self):
        self.write_data()

