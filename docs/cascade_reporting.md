# Is this a waterfall?

No we are not talking waterfalls. We are talking about customising reporting for DST for _M. tuberculosis_ in a clinical and/or PH setting. In these settings, you may want to capture the variants and interpretations of all drugs in a sequence, but due to poor performance of inferring resistance for some drugs or the particulars of the jurisdiction or location for which the report is being generated, it may not be appropriate to report all of these to a client. So `tbtAMR` allows you to generate a "cascade report". These reports are based on rules supplied by the user in the catalogue config file (see [here](https://github.com/MDU-PHL/tbtamr/blob/master/example_criteria/db_config_who_v2.json) for an example and [here](customisation.md)) for explanation. What this means is that the user can specify: 

**If resistance is observed to these drugs, then report these drugs. Otherwise just report default drugs**. 

This is not a feature that would be suited to all, but may be useful in some situations. Note that `tbtAMR` will generate and keep mutliple reports, depending on the setup in your configuration files.