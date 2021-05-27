# Annotation Summary

Version: 0.0.1

## What to use it for?

The Annotation Summary can be used to give you an excellent overview over your data. It will compare your DIAMOND-Megan6-Pipeline annotated metagenomic data to annotated metagenomic data from [MG-Rast](https://www.mg-rast.org/index.html?stay=1).

## Prerequisites

This program is written Python and requires Python 3.6+. Furthermore, the following packages are required:

```
setuptools~=49.6.0
pandas~=1.2.4
matplotlib~=3.4.1
numpy~=1.19.2
scipy~=1.6.2
```

## Options

Several Options are available. For some of them a stable internet connection is required.

### MGM_ID

```
-id MGM_ID, --mgm_id MGM_ID
```

This option can be used when no MG-Rast ID is locally available. If included, it will download the necessary annotation files from the website, directly. The ID should be same as the one that was imported into Megan (The names must not be the same but the metagenomes should be). *Note*: This step might take a while due to very slow downloading speed of the Website.

Also, the `-id` option has a higher priority than `-c` or `-k` which are explained in more detail further below.

### Output

```
-o OUTPUT, --output OUTPUT
```

The Path and file name to save the output to. Output will be a summary csv file comparing all annotations with each other. Also, graphics will be saved to the output path.

### MEGAN_COG

```
-m_c MEGAN_COG, --megan_cog MEGAN_COG
```

Use this to import your exported [Megan6](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/) file. Export the required file from Megan6 by opening the EGGNOG viewer in the main panel. From there, uncollapse the annotation tree (Top-Panel => Tree => Uncollapse All) and select all nodes and leaves (Ctrl/Cmd + a). From the Top-Panel, choose *File* and move the cursor to *Export* and click on *Text (CSV) Format*. Export *eggnogName_to_count* and separator *comma*. The file that is then exported can be used for the analysis.

### MEGAN_KEGG

```
-m_k MEGAN_KEGG, --megan_kegg MEGAN_KEGG
```

Use this to import your exported [Megan6](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/) file. Export the required file from Megan6 by opening the KEGG viewer in the main panel. From there, uncollapse the annotation tree (Top-Panel => Tree => Uncollapse All) and select all nodes and leaves (Ctrl/Cmd + a). From the Top-Panel, choose *File* and move the cursor to *Export* and click on *Text (CSV) Format*. Export *eggnogName_to_count* and separator *comma*. The file that is then exported can be used for the analysis.

**NOTE:** It is possible to analyze both, COG and KEGG in the same run.

### MG-Rast COG

```
-c COG [COG ...], --cog COG [COG ...]
```

Include this option with all functional COG files that can be downloaded from MG-Rast (for a metagenome). If not locally available, provide the MG-Rast ID (mgm...) with this flag to download it similar to the `-id` option. 

**Note:** If both, `-id` and `-c` are included, then `-id` is treated with a higher priority.

### MG-Rast KEGG

```
-k KEGG [KEGG ...], --kegg KEGG [KEGG ...]
```

Include this option with all functional KEGG files that can be downloaded from MG-Rast (for a metagenome). If not locally available, provide the MG-Rast ID (mgm...) with this flag to download it similar to the `-id` option. 

**Note:** If both, `-id` and `-k` are included, then `-id` is treated with a higher priority.

**NOTE:** It is possible to analyze both, COG and KEGG in the same run.

### Temporary Files

#### TEMP_DIR

```
-t TEMP_DIR, --temp_dir TEMP_DIR
```

Use this option to create your own temporary file directory.

#### keep_temp

```
--keep_temp
```

Include this option to keep temporary files. It is recommended to include this option due to the very slow downloading speed of MG-Rast Servers.






