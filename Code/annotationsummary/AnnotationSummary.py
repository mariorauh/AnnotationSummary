# Copyright: Mario Rauh
# Version: 0.0.2
import argparse as ap
import bin.KEGG as KEGG, bin.COG as COG


def command_line():
    parser = ap.ArgumentParser("Metagenomic Data Collection (via MG-Rast) - KOG ID Collector")

    parser.add_argument("-id", "--mgm_id", help="MG-Rast Id", type=str)
    parser.add_argument("-g_c", "--groups_cog", help="Input CSV file containing the COG counts. "
                                                     "Note: The groups need to "
                                                     "be in the same order as the input files. Optional.", type=str)
    parser.add_argument("-g_k", "--groups_kegg", help="Input CSV file containing the KEGG counts. "
                                                      "Note: The groups need to "
                                               "be in the same order as the input files. Optional.", type=str)
    parser.add_argument("-o", "--output", help="The output file's name", default="annotation_summary", type=str)
    parser.add_argument("-m_c", "--megan_cog", help="Input megan export file from EGGNOG/COG:"
                                                    " Name_to_Count. Use in comb. with \"-c\". Optional.", type=str)
    parser.add_argument("-m_k", "--megan_kegg", help="Input megan export file from KEGG:"
                                                     " Name_to_Count. Use in comb. with \"-k\". Optional.", type=str)
    parser.add_argument("-t", "--temp_dir", help="Set a temporary directory to save files in. Will be deleted "
                                                 "afterwards", type=str, default="temp/")
    parser.add_argument("--keep_temp", help="Include to keep temporary files.", action='store_true')
    parser.add_argument("--include_groups", help="Include in combination with MG-Rast ID to download Functional "
                                                 "Groups to include in summary."
                        , type=str)
    parser.add_argument("-c", "--cog", help="Include in combination with cog files to analyze those. If "
                                            "included without files: Add \"mgm...\" "
                                            "=> Adds necessary MG-Rast files based on given id.",
                        nargs='+')
    parser.add_argument("-k", "--kegg", help="Include in combination with kegg files to analyze those. If "
                                            "included without files: Add \"mgm...\" "
                                            "=> Adds necessary MG-Rast files based on given id.",
                        nargs='+')

    return vars(parser.parse_args())


def main():

    print("Annotation Summary")
    args = command_line()
    mgm, temp_dir, cog_files = args["mgm_id"], args["temp_dir"], args["cog"]
    include_groups, output, megan_cog = args["include_groups"], args["output"], \
                                                args["megan_cog"]
    keep_temp, kegg_files, megan_kegg = args["keep_temp"], args["kegg"], args["megan_kegg"]
    groups_cog, groups_kegg = args["groups_cog"], args["groups_kegg"]

    if cog_files is not None:

        print("Analyzing COG Content")
        b = False

        if mgm is not None:
            b, temp_dir = COG.get_functional_cog(mgm, temp_dir)

        if b and len(cog_files) == 1:

            COG.run_with_id(temp_dir, groups_cog, include_groups, mgm, megan_cog, output, keep_temp)

        else:
            COG.run_without_id(temp_dir, groups_cog, include_groups, mgm, megan_cog, output, keep_temp, cog_files)

    temp_dir = args["temp_dir"]

    if kegg_files is not None:
        print("Analyzing KEGG Content")
        b = False

        if mgm is not None:
            b, temp_dir = KEGG.get_functional_kegg(mgm, temp_dir)

        if b and len(kegg_files) == 1:
            KEGG.run_with_id(temp_dir, groups_kegg, include_groups, mgm, megan_kegg, output, keep_temp)

        else:
            KEGG.run_without_id(temp_dir, groups_kegg, include_groups, mgm, megan_kegg, output, keep_temp, kegg_files)


if __name__ == '__main__':
    main()