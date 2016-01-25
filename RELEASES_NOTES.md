# v1.2.0  [dev]
### Update notes:
  The repository structure has been changed: tools have been moved in sub-directories.
  If you update your previous version to this version it is necessary to re-do the following installation steps: 6. Add tools in galaxy and 7. Set memory and parallelisation settings.

### Bugs fixes:
  Fix bug with heatmap hover in affiliations_stat.

### Functions added:
  * Default BIOM becomes mono-line.
  * Reduce the BIOM file size (identation 4 -> 2 in pretty_print).
  * Change the repository structure (one directory by tool).
  * Use biomTools in normalisation.
  * Number of sampled sequences becomes required in normalisation tool.
  * Add better wrapping for sample name fields in preprocess.
  * Add check on minimum amplicon size value in preprocess.


# v1.1.0  [2015-11-30]
### Bugs fixes:
  * Fix bug with empty blast_taxonomy in BIOM.
  * Fix bug with checkall in chrome.
  * Fix bug with CSV extraction in HTML outputs.
  * Fix bug 'Exception in thread QueueFeederThread' in slow systems.

### Functions added:
  * Add '% kept' in preprocess.
  * Add checks on samples names in preprocess.
  * Add DenseData management in BIOM.
  * Add the biom1 datatype in tools.


# v1.0.0  [2015-09-18]
  First package.