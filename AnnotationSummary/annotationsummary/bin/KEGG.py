#file -- KEGG.py --

import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import os
import json

def get_functional_kegg(mgm:str, temp_dir:str):

    # boolean Value b that is set to True if all files were downloaded correctly.
    b = False

    if temp_dir[-1] != "/":
        temp_dir+="/"

    exists = os.path.exists(temp_dir)

    if not exists:
        os.system(f"mkdir {temp_dir}")

    temp_dir+="kegg/"

    exists = os.path.exists(temp_dir)
    if not exists:
        os.system(f"mkdir {temp_dir}")

    try:
        curl1 = f"curl -sS \"https://api-ui.mg-rast.org/annotation/sequence/{mgm}.3?source=KO&" \
                f"type=ontology&filter_level=level1&browser=1&filter=Metabolism\" > {temp_dir}metabolism.csv"
        print(f"Trying to download Metabolism of {mgm} to {temp_dir}metabolism.csv ... ")
        #os.system(curl1)
        #print(f"Successfully downloaded METABOLISM.")
        print(f"Trying to download Genetic Information Processing of {mgm} to {temp_dir}genetics.csv ...")
        curl2 = f"curl -sS \"https://api-ui.mg-rast.org/annotation/sequence/{mgm}.3?source=KO&" \
                f"type=ontology&filter_level=level1&browser=1&filter=Genetic%20Information%20Processing\" > " \
                f"{temp_dir}genetics.csv"
        #os.system(curl)
        #print(f"Successfully downloaded INFORMATION STORAGE AND PROCESSING.")
        print(f"Trying to download Environmental Information Processing of {mgm} to {temp_dir}environment.csv")

        curl3 = f"curl -sS \"https://api-ui.mg-rast.org/annotation/sequence/{mgm}.3?source=KO&" \
                f"type=ontology&filter_level=level1&browser=1&filter=Environmental%20Information%20Processing\" > " \
                f"{temp_dir}environment.csv"
        #os.system(curl)
        #print(f"Successfully downloaded CELLULAR PROCESSES AND SIGNALING.")

        print(f"Trying to download Cellular Processes of {mgm} to {temp_dir}cellular.csv")
        curl4 = f"curl -sS \"https://api-ui.mg-rast.org/annotation/sequence/{mgm}.3?source=KO&type=ontology&" \
                f"filter_level=level1&browser=1&filter=Cellular%20Processes\" > {temp_dir}cellular.csv"
        #os.system(curl4)
        #print(f"Successfully downloaded POORLY CHARACTERIZED.")

        print(f"Trying to download Organismal Systems of {mgm} to {temp_dir}systems.csv")
        curl5 = f"curl -sS \"https://api-ui.mg-rast.org/annotation/sequence/{mgm}.3?source=KO&" \
                f"type=ontology&filter_level=level1&browser=1&filter=Organismal%20Systems\" > {temp_dir}systems.csv"

        print(f"Trying to download Human Diseases of {mgm} to {temp_dir}diseases.csv")
        curl6 = f"curl -sS \"https://api-ui.mg-rast.org/annotation/sequence/{mgm}.3?source=KO&" \
                f"type=ontology&filter_level=level1&browser=1&filter=Human%20Diseases\" > {temp_dir}diseases.csv"

        os.system(f"{curl1} & {curl2} & {curl3} & {curl4} & {curl5} & {curl6} wait")
        b = True

    except:
        print("An error occurred while downloading the necessary mgm files. Please check your input and try again.")
    if b:
        return True, temp_dir
    return b, temp_dir


def to_dataframe(files:list, groups, temp_dir:str):

    print("Merging KEGG files ... ")
    print(f"Checking if functional groups file was given ...")
    g = {}
    if groups is not None:
        print(f"Groups file found.")
        ## to do!
        with open(groups, mode='r') as infile:
            reader = csv.reader(infile)
            keys_values = []
            for row in reader:
                for i in row:
                    keys_values.append(i.split("\t"))


        for i in range(len(keys_values[0])):
            keys = keys_values[0]
            values = keys_values[1]
            g[keys[i]] = values[i]

    # if no mgm id is given:
    if len(files) == 1:
        mgm = files[0]
        b, temp_dir = get_functional_kegg(mgm, temp_dir)
        files = os.listdir(temp_dir)
        filtered_files = []
        for i, file in enumerate(files):
            if ".DS_Store" == file:
                continue
            else:
                filtered_files.append(f"{temp_dir + file}")

        files = filtered_files

    counter = 0
    all_df = []
    for file in files:
        try:
            print(file)
            cols = ["sequence id", "m5nr id (md5sum)", "dna sequence", "semicolon separated list of annotations"]
            df = pd.read_table(file, delimiter="\t", skipinitialspace=True, usecols=cols)
            cols = ["sequence id", "m5nr id (md5sum)", "dna sequence", "semicolon separated list of annotations"]
            #df = pd.read_csv(file, delimiter="\t", usecols=cols)
            df = df.iloc[:-1]
            #print(g.keys())
            if not g:
                temp = ["Not assigned"] * df.shape[0]
            else:
                e = list(g.keys())[counter]
                v = df.shape[0]
                temp = [e] * int(v)

            counter+=1

            #print(len(temp))
            df.insert(loc=4, column="Functional Group", value=temp)

            all_df.append(df)
        except Exception as e:
            print(e)
            print(f"{file} does not exist or can not be loaded correctly.")

    df = pd.concat(all_df, ignore_index=True)
    #print(df.shape[0])
    print(f"Merging completed.")
    return df


def get_ids(df:pd.DataFrame):

    print(f"Checking KEGG IDs ...")
    id_str = df["semicolon separated list of annotations"]
    funcs = df["Functional Group"]
    set_ids = []
    dict_ids = {}
    func_ids = {}
    func_counts = {}
    func_names = {}
    for x,i in enumerate(id_str):
        kegg_id = ""
        b = False

        if i is not np.nan:

            for j in range(11, len(i)):
                if i[j] == "]":
                    break
                else:
                    kegg_id+=i[j]
            func = ""
            for j in range(29, len(i)):
                if i[j] == "]" or i[j] == "[":
                    break
                else:
                    func+=i[j]

            set_ids.append(kegg_id)
            if kegg_id in dict_ids.keys():
                dict_ids[kegg_id] = dict_ids[kegg_id] + 1
            else:
                func_ids[kegg_id] = funcs[x]
                dict_ids[kegg_id] = 1
                func_names[kegg_id] = func

    print(f"Counting Functional Groups ... ")
    for k in funcs:
        if k in func_counts.keys():
            func_counts[k] = func_counts[k] + 1
        else:
            func_counts[k] = 1

    #print(set_ids)
    print(f"Sorted KEGG IDs.")
    return dict_ids, func_ids, func_counts, func_names


def create_csv(count_ids:dict, func_ids:dict, func_counts:dict, func_names:dict, output:str, megan_values):

    print(f"Writing csv to {output}.kegg.csv ... ")
    with open(f"{output}.kegg.csv", "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(["KEGG_ID", "MG_Appearances", "MeganAppearances",
                             "Name"])
        for key in count_ids.keys():
            temp = []
            temp.append(key)
            temp.append(count_ids[key])

            if key in megan_values.keys():
                temp.append(megan_values[key])
            else:
                temp.append(0)

            x = func_ids[key]

            temp.append(func_names[key])
            csv_writer.writerow(temp)

        if megan_values is not None:
            for key in megan_values.keys():
                if key not in count_ids.keys():
                    temp.append(key)
                    temp.append(0)
                    temp.append(megan_values[key])
                    temp.append("Unknown")   # assign unknown functional group because its not part of import csv
                    temp.append("-")
                    temp.append("-")


    print(f"Successfully written {output}.kegg.csv")


def megan_conv(megan:str):
    try:
        print(f"Try loading {megan} file into program ...")

        with open(megan, mode='r') as infile:
            reader = csv.reader(infile)
            megan_values = {}
            for row in reader:

                if row[0][:1] == 'k' or row[0][:1] == 'K':
                    id_val = row[0][:6]
                    count = int(row[1])
                    megan_values[id_val] = count

        print(f"Successfully loaded {megan}.")
        return megan_values

    except:
        print(f"{megan} not found or not the right format.")
        return None


def merge_counts(megan_values:dict, count_ids:dict):

    print(f"Merging Counts and finding Overlapping Counts ...")

    merged_counts = {}
    for ids in megan_values.keys():
        temp = []
        temp.append(megan_values[ids])
        if ids in count_ids.keys():
            temp.append(count_ids[ids])

        merged_counts[ids] = temp

    for ids in count_ids.keys():
        if ids not in merged_counts.keys():
            merged_counts[ids] = [count_ids[ids]]

    #print(merged_counts)
    #print(merged_counts)
    ids_list = []
    mgr_counts = []
    megan_counts = []
    for ids in merged_counts.keys():
        if len(merged_counts[ids]) > 1:
            megan_counts.append(merged_counts[ids][0])
            mgr_counts.append(merged_counts[ids][1])
            ids_list.append(ids)

    print(f"Merged Counts and identified counts.")

    return megan_counts, mgr_counts, ids_list


def create_stacked_bar(megan_counts:list, mgr_counts:list, ids_list:list, output:str):

    print(f"Creating Stacked Bar Chart ... ")

    width = 1.5
    fig,ax = plt.subplots()

    ax.bar(ids_list, megan_counts, width, label='Megan Counts')
    ax.bar(ids_list, mgr_counts, width, bottom=megan_counts, label='MG-Rast Counts')

    ax.set_ylabel('Counts')
    ax.set_xlabel('KEGG IDs')
    ax.legend()

    plt.savefig(f"{output}_stacked.kegg.pdf")
    plt.close()
    print(f"Successfully created Stacked Bar Chart.")

def create_boxplot(megan_counts:list, mgr_counts:list, ids_list:list, output:str):

    print(f"Creating Boxplot ... ")
    all_data = []
    labels=["Megan", "MG-Rast"]
    megan_counts = np.array(megan_counts)
    mgr_counts = np.array(mgr_counts)
    all_data.append(megan_counts)
    all_data.append(mgr_counts)
    plt.boxplot(all_data, vert=True, labels=labels)
    plt.ylabel('Counts')
    plt.title('Boxplot Comparison between Megan Annotation & MG-Rast Annotation')
    plt.savefig(f"{output}_boxplot.kegg.pdf")
    plt.close()
    print(f"Successfully created Boxplot.")


def create_scatter(megan_counts:list, mgr_counts:list, ids_list:list, output:str):

    print(f"Creating Scatter Plot ... ")

    megan_counts = np.array(megan_counts)
    mgr_counts = np.array(mgr_counts)

    plt.plot(megan_counts, mgr_counts,'bo')

    m,b = np.polyfit(megan_counts, mgr_counts, 1)

    plt.plot(megan_counts, m*megan_counts + b)
    plt.ylabel("MG-Rast Counts")
    plt.xlabel("Megan Counts")
    plt.title("Scatter Plot")
    #plt.figtext(.5,.9,'Foo Bar')
    #plt.xticks([0,1000,2000,3000,4000,5000,6000])
    pear = list(pearsonr(megan_counts, mgr_counts))
    plt.figtext(x=0.5, y=0, s=f"r = {round(pear[0],3)}")
    plt.tight_layout()
    plt.savefig(f"{output}_scatter.kegg.pdf")
    plt.close()
    print(f"Successfully created Scatter Plot.")


def remove_temp(temp_dir:str):

    print(f"Removing temporary files from {temp_dir} ... ")
    os.system(f"rm -r {temp_dir}")

def get_groups(mgm:str, temp_dir:str):

    if temp_dir[-1] != "/":
        temp_dir+="/"

    exists = os.path.exists(f"{temp_dir}groups")

    if not exists:
        os.system(f"mkdir {temp_dir}groups")

    print(f"Trying to download functional group file ...")
    try:  #
        curl = f"curl -sS \"https://api-ui.mg-rast.org/metagenome/{mgm}.3?verbosity=stats&detail=ontology\"" \
               f" > {temp_dir}groups/groups.txt"
        # https://api-ui.mg-rast.org/metagenome/mgm4512586.3?verbosity=stats&detail=ontology

        os.system(f"{curl}")
        j = f"{temp_dir}groups/groups.txt"

        print("Successfully downloaded functional group file.")
        data = []
        functions = []

        with open(j, "r") as f:
            data = f.read()

        js = json.loads(data)

        kegg = js["KO"]
        #print(kegg)
        names = []
        num = []
        for k in kegg:
            names.append(k[0])
            num.append(k[1])

        with open(f"{temp_dir}groups/kegg.groups.csv", "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile, delimiter='\t')
            csv_writer.writerow(names)
            csv_writer.writerow(num)

        return f"{temp_dir}groups/kegg.groups.csv"
    except Exception as e:
        print(e)
        return


def run_without_id(temp_dir:str, groups:str, include_groups, mgm:str,
                   megan:str, output:str, keep_temp:bool, files:list):
    print(f"Temporary mgm files given.")

    if groups is None:
        if include_groups is not None:
            groups = get_groups(include_groups, temp_dir)

    df = to_dataframe(files, groups, temp_dir)
    count_ids, func_ids, func_counts, func_names = get_ids(df)

    megan_values = None
    if megan is not None:
        megan_values = megan_conv(megan)
        create_csv(count_ids, func_ids, func_counts, func_names, output, megan_values)
        megan_counts, mgr_counts, ids_list = merge_counts(megan_values, count_ids)
        create_stacked_bar(megan_counts, mgr_counts, ids_list, output)
        create_boxplot(megan_counts, mgr_counts, ids_list, output)
        create_scatter(megan_counts, mgr_counts, ids_list, output)
        # remove temp files
        if not keep_temp:
            remove_temp(temp_dir)

        print(f"\n-------------------- \n")

        print(f"Processed {len(list(megan_values.keys()))} KEGG Functions from Megan")
        print(f"Processed {len(list(count_ids.keys()))} KEGG Functions from MG-Rast\n")

    else:
        create_csv(count_ids, func_ids, func_counts, func_names, output, megan)
        # remove temp files
        if not keep_temp:
            remove_temp(temp_dir)

        print(f"\n-------------------- \n")
        print(f"Processed {len(list(count_ids.keys()))} KEGG Functions from MG-Rast\n")


def run_with_id(temp_dir:str, groups:str, include_groups, mgm:str, megan:str, output:str, keep_temp:bool):

    files = os.listdir(f"{temp_dir}")
    filtered_files = []
    for i, file in enumerate(files):
        if ".DS_Store" == file:
            continue
        else:
            filtered_files.append(f"{temp_dir + file}")
    files = filtered_files

    if groups is None:
        if include_groups is not None:
            groups = get_groups(include_groups, temp_dir)
    df = to_dataframe(files, groups, temp_dir)
    count_ids, func_ids, func_counts, func_names = get_ids(df)


    megan_values = None
    if megan is not None:
        megan_values = megan_conv(megan)
        create_csv(count_ids, func_ids, func_counts, func_names, output, megan_values)
        megan_counts, mgr_counts, ids_list = merge_counts(megan_values, count_ids)
        create_stacked_bar(megan_counts, mgr_counts, ids_list, output)
        create_boxplot(megan_counts, mgr_counts, ids_list, output)
        create_scatter(megan_counts, mgr_counts, ids_list, output)

        # remove temp files
        if not keep_temp:
            remove_temp(temp_dir)

        print(f"\n-------------------- \n")

        print(f"Processed {len(list(megan_values.keys()))} KEGG Functions from Megan")
        print(f"Processed {len(list(count_ids.keys()))} KEGG Functions from MG-Rast\n")

    else:
        create_csv(count_ids, func_ids, func_counts, func_names, output, megan)
        # Remove temporary files
        if not keep_temp:
            remove_temp(temp_dir)

        print(f"\n-------------------- \n")
        print(f"Processed {len(list(count_ids.keys()))} KEGG Functions from MG-Rast\n")




