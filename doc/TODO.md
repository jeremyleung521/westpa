## TODO List for Sphinx Page
* Buttons are different colors
    * this may be adjusted via the html_sidebars option in conf.py
```
# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}
```

* Dates and version numbers
* Descriptive titles
* green/yellow bg color for code (#F6F6F6)
* take out defunct/retired w_tools
* Weighted Ensemble Algorithm section
* WEST>Setup>Running transfer pages to sphinx
* WEST Extension > Plugins MAB scheme? What for?
* Overview is Empty? Delete? Or fill out intro to WESTPA dev
* Update Style Guide - flak8 and black + precommit
* Source Code Management - delete
* WESTPA Modules API? -> python API? (wait to change) 
* tab groups can be changed in /westpa/doc/sphinx_index.rst
* Link back to github wiki, westpa.github.io/westpa
* copyright names and years
* fix faq link on readme landing page, check other links
* side bar command_line_tool_index button should link to page instead of expanding, this might be more navigable


## Regarding autodoc
* westext is not being surveyed by correctly, due to import errors
* westpatools (westpa.cli.tools package) have odd description formatting. Might have to go back to fix docstrings in files.
* Reorganize autodoc structure because all information are 5 levels in



## Completed on 26Jan2021
* Added sidebar captions for (users/developers)
* added checklist and faq pages to users sidebar
* added autosummary sphinx extension
    * autosummary doesn't work too well since most documentation is not appended to each respective function/class
* Organized toc (CLI tools) tree structure to match wiki formatting
* all cli tools have at least some documentation now:
    * w_trace, ploterr, w_crawl, w_fork, w_kinavg, w_kinetics, w_ntop, w_select, w_stateprobs, w_states, and w_succ were added/filled
    * still need to update formatting of these to match previous structure
