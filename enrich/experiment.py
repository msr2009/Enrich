from __future__ import print_function
from enrich_error import EnrichError
from datacontainer import DataContainer
import selection
import time
import sys


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


class Experiment(DataContainer):
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
    results will be based on the ratio of the input and the last timepoint 
    instead of scores.
    """
    def __init__(self, config):
        DataContainer.__init__(self, config)
        self.conditions = dict()
        self.control = None
        self.use_scores = True
        self.normalize_wt = False

        try:
            if 'normalize wt' in config:
                if config['normalize wt'] is True:
                    self.normalize_wt = True
            for cnd in config['conditions']:
                if not cnd['label'].isalnum():
                    raise EnrichError("Alphanumeric label required for condition '%s'" % cnd['label'], self.name)
                for sel_config in cnd['selections']: # assign output base if not present
                    if 'output directory' not in sel_config:
                        sel_config['output directory'] = self.output_base
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

        # ensure consistency for score usage
        if not all(x.use_scores for x in all_selections):
            self.use_scores = False

        # ensure consistency for wild type normalization
        for sel in all_selections:
            sel.normalize_wt = self.normalize_wt



    def calculate(self):
        """
        Calculate scores for all :py:class:`~selection.Selection` objects.
        """
        first = True
        for c in self.conditions:
            s_id = 1
            for s in self.conditions[c]:
                s_label = "%s.%d" % (c, s_id)
                s_id += 1
                s.calculate()
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
                            self.df_dict[dtype] = s.df_dict[dtype][['ratio.%d' % s.timepoints[-1]]]
                            cnames = ["ratio.%s" % s_label]
                        first = False
                    else:
                        for dtype in self.df_dict:
                            self.df_dict[dtype] = self.df_dict[dtype].join(s.df_dict[dtype][['ratio.%d' % s.timepoints[-1]]],
                                how="outer", rsuffix="%s" % s_label)
                            cnames.append("ratio.%s" % s_label)
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
        self.write_data(os.path.join(self.output_base, "experiment_prefilter"))
        # for each filter that's specified
        # apply the filter
        self.filter_stats['total'] = sum(self.filter_stats.values())


    def write_all(self):
        self.write_data()
        for c in self.conditions:
            for sel in self.conditions[c]:
                sel.write_all()



