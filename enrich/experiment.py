from __future__ import print_function
from enrich_error import EnrichError
import selection
import time


def condition_cv_apply_fn(row, condition, use_scores):
    """
    :py:meth:`pandas.DataFrame.apply` function for calculating the 
    coefficient of variation for a variant's score (or ratio) in the 
    condition.
    """
    bc_scores = barcode_data.ix[mapping.variants[row.name]]['score']
    bc_scores = bc_scores[np.invert(np.isnan(bc_scores))]
    cv = stats.variation(bc_scores)
    return pd.Series({'scored.unique.barcodes' : len(bc_scores), \
                      'barcode.cv' : cv})


class Experiment(object):
    """
    Class for a coordinating multiple :py:class:`~.selection.Selection` 
    objects. Creating an 
    :py:class:`~experiment.Experiment` requires a valid *config* object, 
    usually from a ``.json`` configuration file.

    Example config file for a :py:class:`~experiment.Experiment`:

    .. literalinclude:: config_examples/experiment.json


    :download:`Download this JSON file <config_examples/experiment.json>`

    Each experimental ``condition`` contains one or more 
    :py:class:`~.selection.Selection` elements. Multiple 
    :py:class:`~.selection.Selection` objects assigned to the same 
    condition are treated as independent replicates. One ``condition``
    may be specified as a ``control``, and used for control-based score
    corrections.

    .. todo:: Control-based corrections are not yet implemented.

    If any :py:class:`~.selection.Selection` has only two timepoints, the 
    results will be based on ratios of input and selected timepoints 
    instead of scores.
    """
    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.verbose = False
        self.output_directory = None
        self.log = None
        self.conditions = dict()
        self.control = None
        self.df_dict = dict()
        self.use_scores = True

        try:
            self.name = config['name']
            self.output_directory = config['output directory']
            for cnd in config['conditions']:
                if not cnd['label'].isalnum():
                    raise EnrichError("Alphanumeric label required for condition '%s'" % cnd['label'], self.name)
                self.conditions[cnd['label']] = [selection.Selection(x) for x in cnd['selections']]
                if cnd['control']:
                    if self.control is None:
                        self.control = self.conditions[cnd['label']]
                    else:
                        raise EnrichError("Multiple control conditions", self.name)
        except KeyError as key:
            raise EnrichError("Missing required config value %s" % key, 
                              self.name)

        all_selections = list()
        for key in self.conditions:
            all_selections.extend(self.conditions[key])
        for dtype in all_selections[0].df_dict:
            if all(dtype in x.df_dict for x in all_selections):
                self.df_dict[dtype] = True
        if len(self.df_dict.keys()) == 0:
            raise EnrichError("No enrichment data present across all selections", 
                              self.name)

        for key in self.conditions:
            if any(len(x.timepoints) == 2 for x in self.conditions[key]):
                self.use_scores = False


    def enable_logging(self, log):
        """
        Turns on log output for this object. Messages will be sent to the 
        open file handle *log*. 

        .. note:: One log file is usually shared by all objects in the \
        analysis.
        """
        self.verbose = True
        self.log = log
        try:
            print("# Logging started for '%s': %s" % 
                        (self.name, time.asctime()), file=self.log)
        except (IOError, AttributeError):
            raise EnrichException("Could not write to log file", self.name)


    def set_filters(self, config_filters, default_filters):
        """
        Sets the filtering options using the values from the 
        *config_filters* dictionary and *default_filters* dictionary. Filtering 
        options include consistency of scores across replicates.

        .. note:: To help prevent user error, *config_filters* must be a \
        subset of *default_filters*.
        """
        pass


    def calc_selection_scores(self):
        """
        Calculate scores for all :py:class:`~selection.Selection` objects.
        """
        first = True
        for c in self.conditions:
            s_id = 1
            for s in self.conditions[c]:
                s_label = "%s.%d" % (c, s_id)
                s_id += 1
                print("Counting variants for selection %s" % s.name)
                s.calc_all()
                if self.use_scores: # keep the score and r_sq columns
                    if first:
                        for dtype in self.df_dict:
                            self.df_dict[dtype] = s.df_dict[dtype][['score', 'r_sq']]
                            cnames = ["%s.%s" % (x, s_label) for x in ['score', 'r_sq']]
                        first = False
                    else:
                        for dtype in self.df_dict:
                            self.df_dict[dtype] = self.df_dict[dtype].join(s.df_dict[dtype][['score', 'r_sq']],
                                how="outer", rsuffix="%s" % s_label)
                            cnames += ["%s.%s" % (x, s_label) for x in ['score', 'r_sq']]
                else:               # only two timepoints, so keep the ratio
                    if first:
                        for dtype in self.df_dict:
                            self.df_dict[dtype] = s.df_dict[dtype][['ratio.%d' % s.timepoints[1]]]
                            cnames = ["ratio.%s" % s_label]
                        first = False
                    else:
                        for dtype in self.df_dict:
                            self.df_dict[dtype] = self.df_dict[dtype].join(s.df_dict[dtype][['ratio.%d' % s.timepoints[1]]],
                                how="outer", rsuffix="%s" % s_label)
                            cnames.append("ratio.%s" % s_label)
                s.save_data(self.output_directory, clear=True)
        for dtype in self.df_dict:
            self.df_dict[dtype].columns = cnames


    def calc_variation(self):
        """
        Calculate the coefficient of variation for each variant's scores or ratios in each condition.
        """
        for dtype in self.df_dict:
            for c in self.conditions:
                if self.use_scores:
                    c_columns = [x.startswith("score.%s") % c for x in self.df_dict[dtype].columns]
                else:
                    c_columns = [x.startswith("ratio.%s") % c for x in self.df_dict[dtype].columns]
                c_values = self.df_dict[dtype][self.df_dict[dtype].columns[c_columns]]                
                self.df_dict[dtype]['%s.cv' % c] = c_values.apply(stats.variation, axis=1)


    def filter_data(self):
        """
        Apply the filtering functions to the data, based on the filter 
        options present in the configuration object. Filtering is performed 
        using the appropriate apply function.
        """
        self.save_data(os.path.join(self.hdf_dir, "experiment_prefilter"), 
                       clear=False)
        # for each filter that's specified
        # apply the filter
        if self.filters['max barcode variation']:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        barcode_varation_filter, axis=1, 
                        args=[self.filters['max barcode variation']])]
            self.filter_stats['max barcode variation'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min count'] > 0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_count_filter, axis=1, 
                        args=[self.filters['min count']])]
            self.filter_stats['min count'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min input count'] > 0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_input_count_filter, axis=1, 
                        args=[self.filters['min count']])]
            self.filter_stats['min input count'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min rsquared'] > 0.0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_rsq_filter, axis=1, 
                        args=[self.filters['min rsquared']])]
            self.filter_stats['min rsquared'] = \
                    nrows - len(self.df_dict['variants'])

        self.filter_stats['total'] = sum(self.filter_stats.values())
