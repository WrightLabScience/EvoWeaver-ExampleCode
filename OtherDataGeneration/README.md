# Other Data Generation

This folder contains scripts used to generate other data used in this study.

**No scripts or data from this folder are required to make any of the figures**.
Rather, these files were used to generate upstream data (e.g., how KEGG gene indices were retrieved).

**Warning:** These scripts are provided for completeness but may be more difficult to run than
the other scripts in this repository. Many of these analyses were performed on distributed compute
systems. If you have questions on any of the scripts in this folder, feel free to email me (Aidan Lakshman, my email is available on GitHub).

The following folders are included:
- CORUM: additional scripts for CORUM data generation
- KEGG: additional scripts for KEGG data generation (Modules, Complexes)
- STRING: additional scripts for mapping STRING against KEGG. Note that `COG.links.detailed.v12.0.txt` is on Zenodo.

Prior to the first revision of this work, some pairs in the KEGG and STRING were duplicated for a variety of reasons.
These were corrected locally, so a script doesn't exist to trim them.