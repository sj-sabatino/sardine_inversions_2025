# %%
# import sys
import os
import pandas as pd
from operator import attrgetter
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy import stats
import pickle
# lets start by defining a window for pi


class Pi_lineage():
    def __init__(self, lineage_name):
        self.lineage = lineage_name
        self.pi_windows = []
        self.overlapping = None
        self.outliers = []

    def make_region(self, number_of_windows):
        while True:
            starting_window = random.randint(1, len(self.nonoverlapping))
            try:
                stop_window = self.nonoverlapping[starting_window+number_of_windows]

            except IndexError:

                continue
            else:
                if stop_window - self.nonoverlapping[starting_window] == number_of_windows:
                    return self.nonoverlapping[starting_window:starting_window+number_of_windows]

    def make_out_regions(self):
        regions = []
        for number_windows in self.number_window_per_outlier:
            random = self.make_region(number_windows)
            # out_windows += [self.average_pi_across_region(random)]
            regions += [random]
        return regions

    def add_outliers(self, outlier):

        c_outlier = Outlier(scaff=outlier.Chromosome,
                            start=outlier.Start,
                            stop=outlier.Stop,
                            gene=outlier.Gene)

        self.outliers.append(c_outlier)

    def add_pi(self, pi_window):

        c_pi = Pi_window(scaff=pi_window.Chromosome,
                         start=pi_window.Start,
                         stop=pi_window.Stop,
                         pi=pi_window.Pi)
        self.pi_windows.append(c_pi)

    def reset_pi(self):
        self.pi_windows = []

    def sort_outlier_windows(self):
        self.outliers.sort(key=attrgetter('scaff_number', 'start'))

    def sort_pi_windows(self):
        self.pi_windows.sort(key=attrgetter('scaff_number', 'start'))

    def average_pi_across_region(self, region):
        return np.mean([self.pi_windows[x].pi for x in region])

    def pi_across_in_regions(self):
        return [self.average_pi_across_region(x) for x in self.windows_per_in_region]

    def pi_across_out_region(self, regions):
        return [self.average_pi_across_region(x) for x in regions]

    def replicate_pi_across_out_regions(self, n):
        replicates = np.array(self.pi_across_out_region())
        for _ in range(n):
            replicates = np.vstack((replicates, self.pi_across_out_region(self.make_out_regions())))
        return np.mean(replicates, axis=0)

    def overlapping_windows(self):

        self.overlapping = set()
        self.nonoverlapping = set(range(0, len(self.pi_windows)))
        self.number_window_per_outlier = []
        self.windows_per_in_region = []
        for outlier in self.outliers:
            count = 0
            this_region = []
            for index, window in enumerate(self.pi_windows):
                if outlier.overlaps(window):
                    this_region += [index]
                    self.overlapping.add(index)
                    count += 1
                    # if self.pi_windows[index].pi is not None:
                    #    count += 1

            self.windows_per_in_region += [this_region]
            self.number_window_per_outlier += [count]

        self.nonoverlapping = list(self.nonoverlapping.symmetric_difference(self.overlapping))
        self.overlapping = list(self.overlapping)

    def __repr__(self) -> str:
        return f"Lineage {self.lineage} with {len(self.outliers)} outliers"


class Window():
    def __init__(self, scaff, start, stop):
        self.scaff = scaff
        self.scaff_number = int(self.scaff.split("_")[1].split('.1')[0])
        self.start = int(start)
        self.stop = int(stop)

    def overlaps(self, other):
        if self.scaff_number == other.scaff_number and ((self.stop >= other.stop >= self.start) or
                                                        (self.stop >= other.start >= self.start) or
                                                        (self.stop <= other.stop and self.start >= other.start)):
            return True
        else:
            return False
        # this should be fine but do need a pair of tests to make sure

    def __repr__(self):
        return f"{self.scaff} {self.start} {self.stop}"


class Outlier(Window):
    def __init__(self, gene, **kwargs):
        super().__init__(**kwargs)
        self.gene = gene


class Pi_window(Window):
    def __init__(self, pi, **kwargs):
        super().__init__(**kwargs)
        # self.pi = float(pi) if pi != "na" else None
        self.pi = float(pi) if pi != "na" else 0


def tests():

    g_w_1 = Outlier(scaff="scaff_0",
                    start=1,
                    stop=1000,
                    gene="gene1")

    g_w_2 = Outlier(scaff="scaff_1",
                    start=5000,
                    stop=10000,
                    gene="gene2")

    g_w_3 = Outlier(scaff="scaffold_1",
                    start=3310859,
                    stop=3313466,
                    gene="gene3")

    w_w_1 = Pi_window(scaff="scaff_0",
                      start=1,
                      stop=1000,
                      pi="w1")

    w_w_2 = Pi_window(scaff="scaff_0",
                      start=500,
                      stop=1500,
                      pi="w2")

    w_w_3 = Pi_window(scaff="scaff_0",
                      start=1000,
                      stop=1500,
                      pi="w3")

    w_w_4 = Pi_window(scaff="scaff_1",
                      start=1,
                      stop=5000,
                      pi="w4")

    w_w_5 = Pi_window(scaff="scaff_1",
                      start=2000,
                      stop=7000,
                      pi="w5")

    w_w_6 = Pi_window(scaff="scaff_1",
                      start=10001,
                      stop=4857348,
                      pi="w6")

    w_w_7 = Pi_window(scaff="scaffold_1",
                      start=3310000,
                      stop=3330000,
                      pi="w7")
    # w_w_1 overlap exactly with g_w_1
    # w_w_2 overlaps with g_w_1
    # w_w_3 overlaps with the last base of g_w_1
    # w_w_4 overlaps with g_w_2 with first base
    # w_w_5 overlaps with g_w_2
    # w_w_6 does not overlap

    gs = [g_w_1, g_w_2, g_w_3]
    ws = [w_w_1, w_w_2, w_w_3, w_w_4, w_w_5, w_w_6, w_w_7]
    for gene in gs:
        for window in ws:
            string = f"Does Gene {gene.gene} overlap with {window.pi}? {gene.overlaps(window)} "
            print(string)


if __name__ == "__main__":
    import time

    now = time.time()

#    lineages = {"lineage_1": ["kil999", "kil99"],
#                "lineage_2": ["itl999", "itl99"],
#                "lineage_3": ["blk999", "blk99"],
#                "lineage_4": ["alo999", "alo99"]}
    lineages = {"lineage_1": ["inv"]}

    folder = "/Users/sabo/Documents/Projects/sardine/25_pi/inversions"

    excel_file = "/Users/sabo/Documents/Projects/sardine/25_pi/inversions/regions.xlsx"

    with open("./results_regions_average_region_10k.csv", "w") as out, open("results_for_excel.csv", "w") as exc:
        out_line = 1

        out.write("Lineage\tthreshold\tPop\trep\tAveragepi_in\tMedianpi_in\tPiVariance_in\tAveragepi_out\tMedianpi_out\tPiVariance_out\tTtest\tMWU\n")
        for lineage in lineages.keys():
            print(lineage)
            print(lineages.keys())
            files = [x for x in os.listdir(f"{folder}/{lineage}/") if x.endswith(".pi")]

            for threshold in lineages[lineage]:
                print(threshold)
                print(lineages[lineage])

                lineage_obj = Pi_lineage(lineage_name=lineage)
                threshold_outliers = pd.read_excel(excel_file,
                                                   sheet_name=threshold, dtype=str, header=1)
                # threshold_outliers = pd.read_excel("./withinLinsOutliers_01.xlsx",
                #                                   sheet_name=threshold, dtype=str, header=1)

                overlapping_windows = None
                non_overlapping_windows = None

                for index, line in threshold_outliers.iterrows():
                    lineage_obj.add_outliers(line)

                lineage_obj.sort_outlier_windows()
                first = True

                a_lwoverlap = []
                a_lnoverlap = []

                l_lwoverlap = []
                l_lnoverlap = []

                for pi_file in files:
                    landlocked = True if "^" in pi_file else False
                    print(pi_file, landlocked)

                    pi_data = pd.read_csv(f"{folder}/{lineage}/{pi_file}",
                                          sep="\t", names=["Chromosome", "Mid", "a", "b", "Pi"])
                    pi_data['Start'] = pi_data['Mid']-10000
                    pi_data['Stop'] = pi_data['Mid']+10000

                    for index, line in pi_data.iterrows():
                        lineage_obj.add_pi(line)

                    lineage_obj.sort_pi_windows()
                    if lineage_obj.overlapping is None:
                        lineage_obj.overlapping_windows()

                    pi_overlaps = lineage_obj.pi_across_in_regions()

                    average_in = np.mean(pi_overlaps)
                    median_in = np.median(pi_overlaps)
                    variance_in = np.var(pi_overlaps)

                    if first:

                        pi_nonoverlaps_regions = lineage_obj.make_out_regions()
                        out.write("\t".join([str(x) for x in pi_nonoverlaps_regions])+"\n")
                        first = False

                    pi_nonoverlaps = lineage_obj.pi_across_out_region(pi_nonoverlaps_regions)
                    average_out = np.mean(pi_nonoverlaps)
                    median_out = np.median(pi_nonoverlaps)
                    variance_out = np.var(pi_nonoverlaps)

                    tt = stats.ttest_rel(pi_overlaps, pi_nonoverlaps)
                    mwu = stats.mannwhitneyu(pi_overlaps, pi_nonoverlaps)
                    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage,
                                                                                        threshold,
                                                                                        pi_file.split("_")[0],
                                                                                        "1",
                                                                                        average_in,
                                                                                        median_in,
                                                                                        variance_in,
                                                                                        average_out,
                                                                                        median_out,
                                                                                        variance_out,
                                                                                        tt.pvalue,
                                                                                        mwu.pvalue
                                                                                        ))
                    plt.boxplot([pi_overlaps, pi_nonoverlaps])
                    plt.ylim(0.16, 0.25)

                    exc.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage, threshold,
                                                                pi_file.split("_")[0], "overlap", len(pi_overlaps),  "\t".join(map(str, pi_overlaps))))
                    exc.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage, threshold,
                                                                pi_file.split("_")[0], "no_overlap", len(pi_nonoverlaps), "\t".join(map(str, pi_nonoverlaps))))
                    plt.savefig(f'{out_line}_{lineage}_{threshold}_1_{pi_file.split("_")[0]}.png')
                    plt.clf()

                    if landlocked:
                        l_lnoverlap += pi_nonoverlaps
                        l_lwoverlap += pi_overlaps
                    else:
                        a_lnoverlap += pi_nonoverlaps
                        a_lwoverlap += pi_overlaps

                    with open(f'{threshold}_{pi_file.split("_")[0]}.pkl', "wb") as pkl:
                        pickle.dump(lineage_obj, pkl)

                    lineage_obj.reset_pi()

                tt_l = stats.ttest_rel(l_lnoverlap, l_lwoverlap)
                mwu_l = stats.mannwhitneyu(l_lnoverlap, l_lwoverlap)

                tt_a = stats.ttest_rel(a_lnoverlap, a_lwoverlap)
                mwu_a = stats.mannwhitneyu(a_lnoverlap, a_lwoverlap)

                out.write("{}\t{}\t{}\t{}\t{}\n".format(lineage, threshold, "landlocked", tt_l, mwu_l))
                out.write("{}\t{}\t{}\t{}\t{}\n".format(lineage, threshold, "anadromous", tt_a, mwu_a))
                exc.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage, threshold,
                                                            "anadromous", "overlap", len(a_lwoverlap),  "\t".join(map(str, a_lwoverlap))))
                exc.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage, threshold,
                                                            "anadromous", "no_overlap", len(a_lnoverlap), "\t".join(map(str, a_lnoverlap))))
                exc.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage, threshold,
                                                            "landlocked", "overlap", len(l_lwoverlap),  "\t".join(map(str, l_lwoverlap))))
                exc.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lineage, threshold,
                                                            "landlocked", "no_overlap", len(l_lnoverlap), "\t".join(map(str, l_lnoverlap))))

                del(lineage_obj)
