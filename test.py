import pandas as pd
import numpy as np
import itertools
from itertools import combinations
from collections import defaultdict
import re
import json
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import math


np.random.seed(0)


ORDERED_AS = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]


def find_real_codon(row, codons, deg_to_as):

    possible = list(codons.loc[codons.AA == row["AS"], "Codon"])
    deg_c = list(row["degenerated_codon"])
    comb_list = [deg_to_as[x] for x in deg_c]
    combinations = list(itertools.product(*comb_list))
    combinations = list(map("".join, combinations))
    real_codon = [x for x in combinations if x in possible]
    return real_codon


def make_result_df(amino_acids, codon_set, codons, deg_to_as):
    result_df = pd.DataFrame()
    all_amino_acids = [j for i in amino_acids for j in i]
    result_df["AS"] = all_amino_acids

    temp_ls = [[x] * len(y) for x, y in zip(codon_set, amino_acids)]
    temp_ls = [j for i in temp_ls for j in i]
    result_df["degenerated_codon"] = temp_ls
    result_df["deg_to_codon"] = result_df.apply(
        lambda x: find_real_codon(x, codons, deg_to_as)[0], axis=1
    )
    result_df = pd.merge(
        result_df,
        codons[["Codon", "Frac."]],
        how="left",
        left_on="deg_to_codon",
        right_on="Codon",
    )
    result_df.drop("deg_to_codon", axis=1, inplace=True)
    result_df["occ.Freq."] = result_df.degenerated_codon.apply(
        lambda x: decode_codon(x, deg_to_as)
    )

    return result_df


def decode_codon(deg_codon, deg_to_as):
    comb_list = [deg_to_as[x] for x in list(deg_codon)]
    combinations = list(itertools.product(*comb_list))
    combinations = list(map("".join, combinations))
    return 1 / len(combinations)


def find_deg_codons(want_AS, codon_excel):

    codons = (
        codon_excel.copy()
    )  # pd.read_excel(codon_excel)  # pd.read_excel("ecoli_codons.xlsx")

    as_translation = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "*": "STOP",
    }

    codons["AA"] = codons["AA"].map(as_translation)

    iupac_codes = {
        "A": ["R", "W", "M", "D", "H", "V", "N", "A"],
        "G": ["R", "S", "K", "B", "D", "V", "N", "G"],
        "C": ["Y", "S", "M", "B", "H", "V", "N", "C"],
        "T": ["Y", "W", "K", "B", "D", "H", "N", "T"],
    }

    all_codon_dict = dict()
    all_codons_individually = dict()

    for c in codons.AA.unique():
        df_c = codons.loc[codons.AA == c]
        temp_list = []
        all_codons_individually[c] = []

        for cod in df_c.Codon:

            list1 = iupac_codes[cod[0]]
            list2 = iupac_codes[cod[1]]
            list3 = iupac_codes[cod[2]]

            all_codon_combinations = []

            for i in list1:
                for j in list2:
                    for k in list3:
                        all_codon_combinations.append(i + j + k)
            temp_list += all_codon_combinations
            all_codons_individually[c].append(all_codon_combinations)

        all_codon_dict[c] = set(temp_list)

    cant_use_codons = []

    want_codons_dict = {x: [] for x in want_AS}

    for k in list(all_codon_dict.keys()):
        if k in want_AS:
            want_codons_dict[k] = list(set(list(all_codon_dict[k])))
        else:
            cant_use_codons += list(all_codon_dict[k])

    cant_use_codons = list(set(cant_use_codons))

    deg_to_as = defaultdict(list)

    for a in iupac_codes.keys():
        for element in iupac_codes[a]:
            deg_to_as[element].append(a)

    aa = []

    for key in list(want_codons_dict.keys()):
        possible = list(codons.loc[codons.AA == key, "Codon"])
        for c in want_codons_dict[key]:
            comb_list = [deg_to_as[x] for x in list(c)]
            combinations_p = list(itertools.product(*comb_list))
            combinations_p = list(map("".join, combinations_p))
            real_codon = [x for x in combinations_p if x in possible]
            if len(real_codon) > 1:
                cant_use_codons += ["".join(list(c))]

    for k in list(want_codons_dict.keys()):
        want_codons_dict[k] = [
            x for x in want_codons_dict[k] if x not in cant_use_codons
        ]

    A = list(want_codons_dict.values())
    A = list(map(set, A))
    U = set.union(*A)

    codon_to_as = defaultdict(list)

    for a in all_codon_dict.keys():
        for element in all_codon_dict[a]:
            codon_to_as[element].append(a)

    deg_to_as = defaultdict(list)

    for a in iupac_codes.keys():
        for element in iupac_codes[a]:
            deg_to_as[element].append(a)

    result = defaultdict(list)
    best_result = []

    for i in range(1, 5):
        # print(i)
        combs = combinations(U, i)
        for c in combs:
            if all(set(c) & l for l in A):
                ls = [codon_to_as[x] for x in set(c)]
                concat_list = [j for i in ls for j in i]
                result[len(c)].append((set(c), ls, len(concat_list)))

                if len(concat_list) == len(want_AS):
                    best_result.append(make_result_df(
                        ls, set(c), codons, deg_to_as))
        if best_result:
            break

    deg_bases = []
    for l in list(iupac_codes.values()):
        for k in l:
            deg_bases += [x for x in k if x not in ["A", "T", "G", "C"]]

    deg_bases = list(set(deg_bases))

    cheapest = 10
    cheapest_df = defaultdict(list)

    for subdf in best_result:
        all_codons = list(subdf["degenerated_codon"].unique())
        all_codons = list(set([x for j in all_codons for x in j]))
        unique_mixes = [x for x in all_codons if x in deg_bases]

        if cheapest >= len(unique_mixes):
            cheapest = len(unique_mixes)
            cheapest_df[cheapest].append(subdf)

    cheapest_df = cheapest_df[cheapest]
    max_freq = 0

    best_df = pd.DataFrame()

    for subdf in cheapest_df:
        if subdf["Frac."].min() > max_freq:
            max_freq = subdf["Frac."].min()
            best_df = subdf

    best_df["factor"] = best_df["occ.Freq."].max() / best_df["occ.Freq."]

    return best_df, cheapest_df


def find_normal_codons(want_AS, codon_excel):

    codons = (
        codon_excel.copy()
    )  # pd.read_excel(codon_excel)  # pd.read_excel("ecoli_codons.xlsx")

    as_translation = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "*": "STOP",
    }

    codons["AA"] = codons["AA"].map(as_translation)

    codon_list = []
    for amino in want_AS:
        codon_list.append(
            codons.loc[codons.AA == amino]
            .sort_values("Frac.", ascending=False)["Codon"]
            .values[0]
        )
    return codon_list


def calculate_melting_temp(seq_in):
    try:
        return mt.Tm_NN(seq_in, nn_table=mt.DNA_NN1, Na=50, dnac1=50, K=0, dnac2=0)
    except:
        return -1


def reverse_primer_basic(seq, site, configs):

    end_idx = site * 3
    start_idx = max(0, end_idx - configs.length_tolerance_range[-1])

    subseq = seq[start_idx:end_idx]
    variable_part = subseq[: configs.tolerance_range]

    g_t_idxes = np.array(
        [
            y
            for x, y in zip(variable_part, np.arange(len(variable_part)))
            if x in ["G", "C"]
        ]
    )

    melting_temps = np.array(
        [calculate_melting_temp(subseq[x:]) for x in g_t_idxes])

    if any(
        (melting_temps > configs.temperature_tolerance_range[0])
        & (melting_temps < configs.temperature_tolerance_range[1])
    ):
        in_range = melting_temps[
            (melting_temps > configs.temperature_tolerance_range[0])
            & (melting_temps < configs.temperature_tolerance_range[1])
        ]
        g_t_idxes_in_range = g_t_idxes[
            (melting_temps > configs.temperature_tolerance_range[0])
            & (melting_temps < configs.temperature_tolerance_range[1])
        ]
        start_idx = g_t_idxes_in_range[in_range.argmin()]

    elif any(
        (melting_temps > configs.temperature_limit[0])
        & (melting_temps < configs.temperature_limit[1])
    ):
        in_range = melting_temps[
            (melting_temps > configs.temperature_limit[0])
            & (melting_temps < configs.temperature_limit[1])
        ]
        g_t_idxes_in_range = g_t_idxes[
            (melting_temps > configs.temperature_limit[0])
            & (melting_temps < configs.temperature_limit[1])
        ]
        start_idx = g_t_idxes_in_range[in_range.argmin()]

    else:
        start_idx = 0
        # print(f"no suiatable REVERSE primer found at site {site + 1}")

    rev_primer = subseq[start_idx:].reverse_complement()

    return rev_primer, calculate_melting_temp(rev_primer), len(rev_primer)


def forward_primer_basic(seq, site, configs):

    start_idx = (site + 2) * 3
    end_idx = start_idx + configs.length_tolerance_range[-1]

    subseq = seq[start_idx:end_idx]
    variable_part = subseq[-configs.tolerance_range:]

    g_t_idxes = np.array(
        [
            y
            for x, y in zip(
                variable_part, np.arange(
                    len(subseq) - len(variable_part), len(subseq))
            )
            if x in ["G", "C"]
        ]
    )
    melting_temps = np.array(
        [calculate_melting_temp(subseq[: x + 1]) for x in g_t_idxes]
    )

    if any(
        (melting_temps > configs.temperature_tolerance_range[0])
        & (melting_temps < configs.temperature_tolerance_range[1])
    ):
        in_range = melting_temps[
            (melting_temps > configs.temperature_tolerance_range[0])
            & (melting_temps < configs.temperature_tolerance_range[1])
        ]
        g_t_idxes_in_range = g_t_idxes[
            (melting_temps > configs.temperature_tolerance_range[0])
            & (melting_temps < configs.temperature_tolerance_range[1])
        ]
        end_idx = g_t_idxes_in_range[in_range.argmin()]

    elif any(
        (melting_temps > configs.temperature_limit[0])
        & (melting_temps < configs.temperature_limit[1])
    ):
        in_range = melting_temps[
            (melting_temps > configs.temperature_limit[0])
            & (melting_temps < configs.temperature_limit[1])
        ]
        g_t_idxes_in_range = g_t_idxes[
            (melting_temps > configs.temperature_limit[0])
            & (melting_temps < configs.temperature_limit[1])
        ]
        end_idx = g_t_idxes_in_range[in_range.argmin()]

    else:
        end_idx = start_idx + configs.length_tolerance_range[-1]
        # print(f"no suiatable FORWARD primer found at site {site + 1}")

    forw_primer = subseq[: end_idx + 1]

    return forw_primer, calculate_melting_temp(forw_primer), len(forw_primer)


def silent_mutation(silent_codon, configs):

    if configs.silent == False:
        return silent_codon
    print("does this run?")

    codon_df = configs.codons_file.copy()  # pd.read_excel(configs.codons_file)

    AS = codon_df.loc[codon_df.Codon == silent_codon, "AA"].values[0]
    index_AS = list(codon_df.loc[codon_df.Codon ==
                    str(silent_codon), "AA"].index)

    alternative_codons = codon_df.loc[codon_df.AA == str(AS)]
    alternative_codons = [
        x for x in alternative_codons.index if x not in index_AS]
    sub_df = codon_df.iloc[alternative_codons].reset_index(drop=True)
    if len(sub_df) >= 1:
        sub_df = sub_df.iloc[sub_df["Frac."].idxmax()].reset_index(drop=True)

    replacement = sub_df[0] if len(sub_df) == 4 else silent_codon
    return replacement


def design_basic_primer(seq, site, degen_codons, configs):

    site = site - 1
    fw, fw_t, fw_l = forward_primer_basic(seq, site, configs)
    rv, rv_t, rv_l = reverse_primer_basic(seq, site, configs)

    silent_site = seq[(site + 1) * 3: (site + 2) * 3]
    silent_site = silent_mutation(silent_site, configs)

    primer_dict = dict()

    primer_dict["reverse" + str(site + 1)] = [rv, rv_t, rv_l]

    for codon in degen_codons:
        t_p = rv.reverse_complement() + codon + silent_site + fw
        primer_dict["forward_" + codon + "_" +
                    str(site + 1)] = [t_p, fw_t, fw_l]

    return primer_dict


def design_quickchange_primer(seq, site, degen_codons, configs):

    site = site - 1
    fw, fw_t, fw_l = forward_primer_basic(seq, site, configs)
    rv, rv_t, rv_l = reverse_primer_basic(seq, site - 1, configs)

    silent_site_f = seq[(site + 1) * 3: (site + 2) * 3]
    silent_site_f = silent_mutation(silent_site_f, configs)

    silent_site_r = seq[(site - 1) * 3: (site) * 3]
    silent_site_r = silent_mutation(silent_site_r, configs)

    primer_dict = dict()
    for codon in degen_codons:
        forward_p = rv.reverse_complement() + silent_site_r + codon + silent_site_f + fw
        reverse_p = forward_p.reverse_complement()
        primer_dict["rv" + "_" + codon + "_" +
                    str(site + 1)] = [reverse_p, rv_t, rv_l]
        primer_dict["fw" + "_" + codon + "_" +
                    str(site + 1)] = [forward_p, fw_t, fw_l]

    return primer_dict


def combinatorial_primers(seq, sites, degen_codons, is_quickchange=True, configs=None):

    sites = np.array(sites) - 1

    if configs.silent == False:
        fw, fw_t, fw_l = forward_primer_basic(seq, sites.max() - 1, configs)
        fw = fw[3:]
    else:
        fw, fw_t, fw_l = forward_primer_basic(seq, sites.max(), configs)

    if not is_quickchange:
        rv, rv_t, rv_l = reverse_primer_basic(seq, sites.min(), configs)

    elif is_quickchange and configs.silent:
        rv, rv_t, rv_l = reverse_primer_basic(seq, sites.min() - 1, configs)

    else:
        rv, rv_t, rv_l = reverse_primer_basic(seq, sites.min(), configs)
        rv = rv[3:]

    silent_sites_list_fv = []

    for s in sites:
        silent_site = seq[(s + 1) * 3: (s + 2) * 3]
        silent_site = silent_mutation(silent_site, configs)
        silent_sites_list_fv.append(silent_site)

    silent_sites_list_rv = []

    if is_quickchange:
        for s in sites:
            silent_site_r = seq[(s - 1) * 3: (s) * 3]
            silent_site_r = silent_mutation(silent_site_r, configs)
            silent_sites_list_rv.append(silent_site_r)

    else:
        for s in sites:
            silent_site_r = seq[(s - 1) * 3: (s) * 3]
            silent_sites_list_rv.append(silent_site_r)

    if is_quickchange:
        forward_p = rv.reverse_complement(
        ) + silent_sites_list_rv[0] + degen_codons[0]
    else:
        forward_p = rv.reverse_complement() + degen_codons[0]

    for i in range(len(sites) - 1):

        g_b = sites[i] + 2
        g_e = sites[i + 1] - 1

        gap = seq[g_b * 3: g_e * 3]

        diff = sites[i + 1] - sites[i]

        if diff == 1:
            forward_p += degen_codons[i + 1]

        elif diff == 2:
            forward_p += silent_sites_list_fv[i] + degen_codons[i + 1]

        elif diff == 3:
            forward_p += (
                silent_sites_list_fv[i]
                + silent_sites_list_rv[i + 1]
                + degen_codons[i + 1]
            )

        else:
            forward_p += silent_sites_list_fv[i]
            forward_p += gap
            forward_p += silent_sites_list_rv[i + 1] + degen_codons[i + 1]

    forward_p += silent_sites_list_fv[-1] + fw

    if is_quickchange:
        reverse_p = forward_p.reverse_complement()
    else:
        reverse_p = rv

    primer_dict = dict()

    mgs_fv = "fv_"
    mgs_rv = "rv_"

    for i in range(len(sites)):
        mgs_fv += str(degen_codons[i]) + "_" + str(sites[i] + 1)
        mgs_rv += str(degen_codons[i]) + "_" + str(sites[i] + 1)

        if i < len(sites) - 1:
            mgs_fv += "_"
            mgs_rv += "_"

    primer_dict[mgs_fv] = [forward_p, fw_t, fw_l]
    primer_dict[mgs_rv] = [reverse_p, rv_t, rv_l]

    return primer_dict


def prepare_dict(site_codon_dict, configs):

    updated_dict = {}
    degen_excels = defaultdict(lambda: 1)

    previousloop = ""
    previous_loop_codons = ""

    for k, v in site_codon_dict.items():
        newval = ""
        v = v.split("/")
        k_split = k.split("/")
        newvalue = []

        for cntr, s in enumerate(v):
            news = s.split(",")
            updated_news = []
            for lastloop in news:
                if "+" in lastloop:
                    lastloop = list(lastloop)[1:]
                    loop_codons = find_deg_codons(
                        lastloop, configs.codons_file)[0]

                    degen_excels.update(
                        dict(loop_codons["degenerated_codon"].value_counts())
                    )

                    loop_codons = list(
                        loop_codons["degenerated_codon"].unique())
                    updated_news.append(loop_codons)
                elif "-" in lastloop:
                    lastloop = list(lastloop)[1:]
                    if lastloop == previousloop:
                        loop_codons = previous_loop_codons
                    else:
                        loop_codons = find_normal_codons(
                            lastloop, configs.codons_file)
                        previousloop = lastloop
                        previous_loop_codons = loop_codons

                    updated_news.append(loop_codons)

                elif "#WT" in lastloop:
                    new_s_ = int(k_split[cntr]) - 1
                    wt_codon = configs.myseq[new_s_ * 3: new_s_ * 3 + 3]
                    updated_news.append([str(wt_codon)])

                elif "#19c" in lastloop:
                    new_s_ = int(k_split[cntr]) - 1
                    wt_codon = configs.myseq[new_s_ * 3: new_s_ * 3 + 3]
                    wt_res = wt_codon.translate()

                    remaining_res = [
                        x
                        for x in ORDERED_AS
                        if x not in wt_res
                    ]
                    loop_codons = find_normal_codons(
                        remaining_res, configs.codons_file)
                    updated_news.append(loop_codons)

                else:
                    updated_news.append([lastloop])

            updated_news = [i_ for b_ in updated_news for i_ in b_]
            updated_news = ",".join(updated_news)
            newvalue.append(updated_news)

        newvalue = "/".join(newvalue)
        updated_dict[k] = newvalue

    return updated_dict, degen_excels


def design_primers(site_codon_dict, configs):

    all_primers = {}

    for k, v in site_codon_dict.items():

        k = list(map(int, k.split("/")))
        v = v.split("/")
        v = [x.split(",") for x in v]
        v = list(itertools.product(*v))

        for i in v:
            i = [x for x in i]
            out_dict = combinatorial_primers(
                configs.myseq, k, i, is_quickchange=configs.quickchange, configs=configs
            )
            all_primers.update(out_dict)

    return all_primers


def get_mixing_ratio(primer_string, degenerate_ratio_dict):

    codons_primer = re.findall("[A-Z][A-Z][A-Z]", primer_string)
    ratio = 1
    for k in codons_primer:
        ratio *= degenerate_ratio_dict[k]

    return ratio


def create_result_df_template(r_dict, fullDf=False):
    results = r_dict.copy()

    primer_df = pd.DataFrame(results).T
    primer_df.rename(
        columns={0: "seq", 1: "melting_temp", 2: "length"}, inplace=True)
    primer_df["seq"] = primer_df["seq"].apply(lambda x: "".join(x))
    print(len(primer_df), "LEN BEFORE")

    if fullDf:
        fv = [x.startswith("fv") for x in primer_df.index]
        rv = [not x for x in fv]
        primer_df = pd.concat(
            [primer_df[fv], primer_df[rv].drop_duplicates("seq")])
        primer_df = primer_df.reset_index().rename(columns={"index": "name"})

    else:
        primer_df = (
            primer_df.drop_duplicates("seq").reset_index().rename(
                columns={"index": "name"})
        )

    print(len(primer_df), "LEN AFTER")
    return primer_df


def create_primer_overlap(primer_df):
    oligo_primers = {}

    for (r, row) in primer_df.iterrows():
        if row["name"].startswith("rv"):
            new_name_base = row["name"].split("_")
            new_name_base = new_name_base[2] + \
                "_" + new_name_base[-1] + "_start"

            _seq = Seq(row["seq"])

            oligo_primers["rv_" + new_name_base] = [
                _seq,
                row["melting_temp"],
                row["length"],
            ]
            oligo_primers["fv_" + new_name_base] = [
                _seq.reverse_complement(),
                row["melting_temp"],
                row["length"],
            ]

        else:
            new_name_base = row["name"].split("_")
            new_name_base = new_name_base[2] + "_" + new_name_base[-1] + "_end"

            _seq = Seq(row["seq"][-row["length"]:])

            oligo_primers["fv_" + new_name_base] = [
                _seq,
                row["melting_temp"],
                row["length"],
            ]
            oligo_primers["rv_" + new_name_base] = [
                _seq.reverse_complement(),
                row["melting_temp"],
                row["length"],
            ]

    return oligo_primers


def postprocess_primers(results_in):
    results = results_in.copy()

    # remove duplicates
    primer_df = create_result_df_template(results, True)
    results_amp = create_primer_overlap(primer_df)

    primer_df = primer_df[
        primer_df.apply(lambda x: not x["name"].startswith("rv"), axis=1)
    ]
    primer_df["name"] = primer_df["name"].apply(
        lambda x: x.replace("fv", "oligo"))
    primer_df["melting_temp"] = "-"
    primer_df["length"] = primer_df["seq"].apply(lambda x: len(x))

    primer_amp_df = create_result_df_template(results_amp)
    primer_amp_df["idx"] = "-"
    primer_amp_df["mut"] = "-"

    primer_df["idx"] = primer_df["name"].apply(lambda x: x.split("_")[4])

    primer_df["mut"] = ORDERED_AS * int(len(primer_df) / 20)
    concat_df = pd.concat([primer_amp_df, primer_df])  # .set_index("name")

    return concat_df


def np_encoder(object):

    if isinstance(object, np.generic):
        return object.item()


class dotdict(dict):
    """dot.notation access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


defaults = {
    "length_tolerance_range": ["13", "26"],
    "temperature_tolerance_range": ["58", "61"],
    "temperature_limit": ["52", "64"],
    "tolerance_range": 15,
    "silent": False,
    "quickchange": False,
    "myseq": "",
    "sites": "",
    "inp": "",
    "codons_file": pd.read_excel("codon_tables/ecoli_codons.xlsx"),
    "outfile": ".",
}


def oligo_design(**args):

    defaults.update(args)
    FLAGS = dotdict(defaults)

    FLAGS.length_tolerance_range = [int(x)
                                    for x in FLAGS.length_tolerance_range]
    FLAGS.temperature_tolerance_range = [
        float(x) for x in FLAGS.temperature_tolerance_range
    ]
    FLAGS.temperature_limit = [float(x) for x in FLAGS.temperature_limit]

    FLAGS.tolerance_range = (
        FLAGS.length_tolerance_range[-1] - FLAGS.length_tolerance_range[0]
    )

    FLAGS.myseq = Seq(FLAGS.myseq.upper())

    site_codon_dict = {
        s: i for (s, i) in zip(FLAGS.sites.split(";"), FLAGS.inp.split(";"))
    }
    site_codon_dict, degen_excels = prepare_dict(site_codon_dict, FLAGS)

    primers = design_primers(site_codon_dict, FLAGS)
    # return primers
    primers = postprocess_primers(primers)

    return primers


def get_anchors(group, grps):
    return grps.loc[group]["min"] - 1, grps.loc[group]["max"] + 1


def create_sites_and_inputs(base_sequence, max_len, buffer=10):
    allAS = "".join(ORDERED_AS)
    protein_sequcene = Seq(base_sequence).translate()

    buffer_beginning_end = buffer

    rn = list(
        range(1 + buffer_beginning_end,
              len(protein_sequcene) - buffer_beginning_end)
    )

    rn_df = pd.DataFrame(rn, columns=["site"])
    rn_df["residues"] = allAS

    num_groups = math.ceil(
        (len(base_sequence) - buffer_beginning_end * 3 * 2) / max_len
    )  # *2 for beginning and end, *3 for base-> codon
    rn_df["group"] = pd.cut(rn_df["site"], num_groups,
                            labels=range(num_groups))

    grps = rn_df.groupby("group")["site"].agg(["min", "max"])

    positions = []
    mutations = []

    for r, row in rn_df.iterrows():

        l_a, u_a = get_anchors(row["group"], grps)

        posi_ = str(l_a) + "/" + str(row["site"]) + "/" + str(u_a)  # + ";"
        muti_ = "#WT/-" + row["residues"] + "/#WT"

        positions.append(posi_)
        mutations.append(muti_)

    mutations = ";".join(mutations)
    positions = ";".join(positions)

    return mutations, positions


def prepare_and_design(bp_sequence, max_length):
    mutations, positions = create_sites_and_inputs(bp_sequence, max_length)
    return oligo_design(myseq=bp_sequence, sites=positions, inp=mutations)


if __name__ == "__main__":
    test_seq = "ATGGCAGAAGCAGCACAGAGCGTTGACCAGCTGATTAAAGCACGTGGTAAAGTTTATTTTGGTGTTGCCACCGATCAGAATCGTCTGACCACCGGTAAAAATGCAGCAATTATTCAGGCAGATTTTGGTATGGTTTGGCCTGAAAATAGCATGAAATGGGATGCAACCGAACCGAGCCAGGGCAATTTTAACTTTGCCGGTGCAGATTATCTGGTTAATTGGGCACAGCAGAATGGTAAACTGATTGGTGGTGGTATGCTGGTTTGGCATAGCCAGCTGCCGAGCTGGGTTAGCAGCATTACCGATAAAAACACCCTGACCAATGTGATGAAAAACCATATCACCACACTGATGACCCGCTATAAAGGTAAAATTCGTGCATGGGATGTTGTGGGTGAAGCCTTTAATGAAGATGGTAGCCTGCGTCAGACCGTTTTTCTGAATGTTATTGGCGAAGATTATATCCCGATTGCATTTCAGACCGCACGTGCAGCAGATCCGAATGCAAAACTGTATATCATGGATTATAACCTGGATAGCGCAAGCTATCCGAAAACACAGGCAATTGTTAATCGTGTTAAACAGTGGCGTGCAGCCGGTGTTCCGATTGATGGTATTGGTAGTCAGACCCATCTGAGCGCAGGTCAAGGTGCGGGTGTTCTGCAGGCACTGCCGCTGCTGGCAAGCGCAGGTACACCGGAAGTTAGCATTCTGATGCTGGATGTTGCAGGCGCAAGCCCGACCGATTATGTTAATGTTGTTAATGCATGCCTGAATGTTCAGAGCTGTGTTGGTATTACCGTTTTTGGTGTGGCAGATCCGGATAGCTGGCGTGCAAGCACCACACCGCTGCTGTTTGATGGTAATTTCAATCCGAAACCGGCATATAATGCCATTGTTCAGGATCTGCAGCAGGGTAGCATTGAAGGTCGTGGTCATCATCACCATCATCATTAA"
    mutations, positions = create_sites_and_inputs(test_seq, 180)

    # mutations, positions = "#WT/-ACDEFGHIKLMNPQRSTVWY/#WT", "265/308/309"

    test_inputs = {
        "myseq": "ATGGCAGAAGCAGCACAGAGCGTTGACCAGCTGATTAAAGCACGTGGTAAAGTTTATTTTGGTGTTGCCACCGATCAGAATCGTCTGACCACCGGTAAAAATGCAGCAATTATTCAGGCAGATTTTGGTATGGTTTGGCCTGAAAATAGCATGAAATGGGATGCAACCGAACCGAGCCAGGGCAATTTTAACTTTGCCGGTGCAGATTATCTGGTTAATTGGGCACAGCAGAATGGTAAACTGATTGGTGGTGGTATGCTGGTTTGGCATAGCCAGCTGCCGAGCTGGGTTAGCAGCATTACCGATAAAAACACCCTGACCAATGTGATGAAAAACCATATCACCACACTGATGACCCGCTATAAAGGTAAAATTCGTGCATGGGATGTTGTGGGTGAAGCCTTTAATGAAGATGGTAGCCTGCGTCAGACCGTTTTTCTGAATGTTATTGGCGAAGATTATATCCCGATTGCATTTCAGACCGCACGTGCAGCAGATCCGAATGCAAAACTGTATATCATGGATTATAACCTGGATAGCGCAAGCTATCCGAAAACACAGGCAATTGTTAATCGTGTTAAACAGTGGCGTGCAGCCGGTGTTCCGATTGATGGTATTGGTAGTCAGACCCATCTGAGCGCAGGTCAAGGTGCGGGTGTTCTGCAGGCACTGCCGCTGCTGGCAAGCGCAGGTACACCGGAAGTTAGCATTCTGATGCTGGATGTTGCAGGCGCAAGCCCGACCGATTATGTTAATGTTGTTAATGCATGCCTGAATGTTCAGAGCTGTGTTGGTATTACCGTTTTTGGTGTGGCAGATCCGGATAGCTGGCGTGCAAGCACCACACCGCTGCTGTTTGATGGTAATTTCAATCCGAAACCGGCATATAATGCCATTGTTCAGGATCTGCAGCAGGGTAGCATTGAAGGTCGTGGTCATCATCACCATCATCATTAA",
        "sites": positions,
        "inp": mutations,
    }

    rr = oligo_design(
        myseq=test_inputs["myseq"],
        sites=test_inputs["sites"],
        inp=test_inputs["inp"],
        codons_file=pd.read_excel(
            "/Users/davidpatsch/programming/oligo/codon_tables/ecoli_codons.xlsx")
    )
    rr.to_csv("test.csv")
