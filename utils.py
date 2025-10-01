# Author: Thomas Cokelaer, 2025
import time
import os
import glob
from tqdm import tqdm
import pandas as pd
from pylab import savefig, tight_layout, legend, xlim, xlabel, ylabel, axvline, xticks, yticks, gca
from sequana import checkm, FastA, FastQ, BUSCO
import subprocess
import wget


genome_size = {
    "veillonella": 2146482,
    "streptococcus" :  2210410,
    "bacteroide": 5205140,
    "cyanobacteria": 4.4e6,
    "leishmania": 33.5e6
}

def download_data(url, dest_dir="data/"):
    os.makedirs(dest_dir, exist_ok=True)
    filename = os.path.join(dest_dir, os.path.basename(url))
    if os.path.exists(filename):
        print(f"✅ {filename} already exists, skipping.")
    else:
        print(f"Downloading {url} to {filename} ...")
        wget.download(url, out=filename)  # shows progress bar
        print(f"✅ Download complete: {filename}")

def download_SRR_with_progress(accession, outdir, expected_size=None):
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{accession}.fastq.gz")

    # Skip if already exists
    if os.path.exists(outfile):
        print(f"✅ {outfile} already exists, skipping.")
        return

    # Launch fastq-dump
    cmd = f"fastq-dump --gzip -O {outdir} {accession}"
    print(f"🚀 Running: {cmd}")
    process = subprocess.Popen(cmd, shell=True)

    # Track progress
    if expected_size:
        pbar = tqdm(total=expected_size, unit="B", unit_scale=True)
        while process.poll() is None:
            if os.path.exists(outfile):
                size = os.path.getsize(outfile)
                pbar.n = size
                pbar.refresh()
            time.sleep(5)
        pbar.close()
    else:
        process.wait()

    if process.returncode == 0:
        print(f"✅ Download complete: {outfile}")
    else:
        print("❌ fastq-dump failed!")



def saveall(filename, dpi=200):
    outdir = os.path.dirname(filename)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Save in all formats
    for ext in ["png", "pdf", "eps"]:
        outfile = f"{filename}.{ext}"
        savefig(outfile, dpi=dpi, bbox_inches="tight")
    print(f"💾 Saved {outfile} in pdf/eps/png")

def to_seconds(x):
    if "-" in x:
        D = x.split("-")[0]
        x = x.split("-")[1]
    else:
        D = 0
    H, M, S = x.split(":")

    return int(S) + int(M)*60 + int(H)*3600 + int(D)*3600*60


class BUSCOFactory:
    def __init__(self):

        self.colors = {
            "Complete": "green",
            "Missing": "red",
            "Duplicated": "orange",
            "Fragmented": "pink"
        }
        self.get_busco_data()

    def get_busco_data(self):

        dfs = []
        paths = list(glob.glob(f"metadata_busco_checkm/*/*/full_table.tsv"))

        for path in tqdm(paths):
            try:
                b = BUSCO(path)
                _, exp, assembler,_ = path.split("/")
                df = b.summary()
                df['exp'] = exp
                if assembler == 'lora':
                    df['assembler'] = 'canu'
                else:
                    df['assembler'] = assembler.replace("_nocirc", "").replace("_circ", " circ")
                df['circlator'] = False if 'nocirc' in assembler else True
                dfs.append(pd.Series(df))
            except Exception as err:
                print(path)
                print(err)
        self.busco_df = pd.concat(dfs ,axis=1).T

    def plot_busco_summary_bar(self, name):
        df_grouped = self.busco_df.query("exp==@name").groupby("assembler")[["S_pc", "D_pc", "F_pc", "M_pc"]].sum()
        df_grouped = df_grouped.iloc[::-1]
        df_grouped = df_grouped.rename(
        columns={"F_pc": "Fragmented", "D_pc": "Duplicated", "S_pc": "Complete", "M_pc": "Missing"}
    )
        df_grouped.plot(
            kind="barh",
            stacked=True,
            figsize=(8, 6),
            fontsize=12,
            color = [self.colors[col] for col in df_grouped.columns if col!= "assembler"],
            width=0.8
        )
        legend(loc="upper left")
        xlim([0,100])
        _ = xlabel("Percentage (%)", fontsize=16)
        _ = ylabel("Assembler", fontsize=16)
        tight_layout()


class CheckMFactory:

    def __init__(self):
        self.colors = {
            "Completeness": "green",
            "Contamination": "orange",
            "Strain heterogeneity": "pink",
        }

    def plot_checkm_summary_bar(self, name):
        assemblers = [x.split("/")[2] for x in glob.glob(f'metadata_busco_checkm/{name}/*//results.txt')]
        #samples = [x.split("/")[3] for x in glob.glob(f'final_analysis/{name}/*/*/checkm/results.txt')]

        c = checkm.MultiCheckM(glob.glob(f'metadata_busco_checkm/{name}/*/results.txt'))
        c.df.columns = [x.replace("_nocirc", "").replace("_circ", " circ") for x in assemblers]

        self.c = c
        df_grouped = c.df.T
        #df_grouped["sample"] = samples
        df_grouped['assembler'] = df_grouped.index
        df_grouped = df_grouped[["assembler", "Completeness", "Contamination" , "Strain heterogeneity"]]
        df_grouped.groupby("assembler")[["Completeness", "Contamination", "Strain heterogeneity"]].sum()
        df_grouped = df_grouped.iloc[::-1]

        df_grouped.plot(
            kind="barh",
            stacked=True,
            figsize=(8, 6),
            fontsize=12,
            color = [self.colors[col] for col in df_grouped.columns if col!= "assembler"],
            width=0.8
        )
        legend(loc="upper left")
        xlim([0,100])
        _ = xlabel("Percentage (%)", fontsize=16)
        _ = ylabel("Assembler", fontsize=16)
        tight_layout()


class PlotContigs:
    def __init__(self, name):
        self.name = name
        filename = f"metadata_contigs/{name}.contigs.csv"
        if os.path.exists(filename) is True:
            print("Reading file")
            self.df = pd.read_csv(filename)
        else:
            print("Creating and saving file. You need FASTA files from final analysis. Not provided here but on Zenodo DOI in the github.com/cokelaer/paper_LORA.")
            self.df = self._get_contigs_df(name)
            self.df.to_csv(filename, index=False, sep=",")

    def _get_contigs_df(self, name):
        assemblers = [x.split("/")[2] for x in glob.glob(f'final_analysis/{name}/*/*/sorted_contigs/*.fasta')]
        samples = [x.split("/")[3] for x in glob.glob(f'final_analysis/{name}/*/*/sorted_contigs/*fasta')]
        filenames = glob.glob(f'final_analysis/{name}/*/*/sorted_contigs/*fasta')

        data = []
        for filename in filenames:
            f = FastA(filename)
            N = sorted(f.get_lengths_as_dict().values())
            sample = filename.split("/")[2]
            assembler = filename.split("/")[3]
            for x in N:
                data.append([sample, assembler, x])

        df = pd.DataFrame(data)
        df.columns = ["assembler", "sample", "length"]

        if name == "leishmania":
            for assembler in ['canu','flye','hifiasm','necat','pecat','unicycler']:
                missing = f"{assembler}" 
                if missing not in df['assembler'].values:
                    df.loc[len(df)] = [missing, "none", 0]
        else:
            for assembler in ['canu','flye','hifiasm','necat','pecat','unicycler']:
                for x in ['circ','nocirc']:
                    missing = f"{assembler}_{x}" 
                    if missing not in df['assembler'].values:
                        df.loc[len(df)] = [missing, "none", 0]
            df['assembler'] = [x.replace("_nocirc", "").replace("_circ", " circ.") for x in df['assembler']]
        return df

    def plot_assembly_results(self):
        
        df_sorted = self.df.sort_values(["assembler", "length"], ascending=[True, False])

        pivot = (
            df_sorted
            .assign(rn=df_sorted.groupby("assembler").cumcount())  # position of each value
            .pivot(index="assembler", columns="rn", values="length")
            .fillna(0)
        )

        # plot stacked bars
        pivot.plot(kind="barh", zorder=10,
                stacked=True, 
                figsize=(8,5),
                legend=False,
                edgecolor="none" )
        xlabel("Total contig length (bp)", fontsize=16)
        axvline(genome_size[self.name.split("_")[0]], ls="--", color="k", zorder=20)
        ylabel("Assembler method", fontsize=16)
        _ = gca().set_yticklabels(gca().get_yticklabels(), ha='left')
        gca().tick_params(axis='y', pad=80)
        _ = xticks(fontsize=12)
        _ = yticks(fontsize=12)
        gca().invert_yaxis()

        ax = gca()
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())  # synchronize y-axis
        ax2.set_yticks(range(len(pivot)))
        ax2.set_yticklabels([(pivot > 0).sum(axis=1)[name] for name in pivot.index])
        ax2.tick_params(axis='y', length=0)  # hide tick marks
        ax2.set_ylabel('Count')

        ax.set_axisbelow(True)     # draw grid under bars
        ax.grid(True, axis="y", linestyle="-", color="0.8")  # style as you like
        tight_layout()

        return df_sorted

def download_veillonella_data():
    os.makedirs("data/veillonella", exist_ok=True)

    if os.path.exists("data/veillonella/m54091_180306_141024.subreads.fastq.gz") is False:
        # You can run fastq-dump ERR3958992 to get the fastq file, or for this notebook simply call:
        download_data("https://zenodo.org/records/13306684/files/m54091_180306_141024.subreads.fastq.gz", "data/veillonella")
    else:
        print("m54091_180306_141024.subreads.fastq.gz already present")

    if os.path.exists("data/veillonella/veillonella.ccs.fastq.gz") is False:
        # for the CCS file, call:
        download_data("https://zenodo.org/records/13306684/files/veillonella.ccs.fastq.gz", "data/veillonella")
    else:
        print("veillonella.ccs.fastq.gz already present")

def download_cyanobacteria_data():
    os.makedirs("data/cyanobacteria", exist_ok=True)

    if os.path.exists("data/cyanobacteria/PCC6711.fastq.gz") is False:
        # You can run fastq-dump ERR3958992 to get the fastq file, or for this notebook simply call:
        download_data("https://zenodo.org/records/17229238/files/PCC6711.fastq.gz", "data/cyanobacteria")


def download_streptococcus_data():
    os.makedirs("data/streptococcus", exist_ok=True)

    if os.path.exists("data/streptococcus/hifi.ccs.fastq.gz") is False:
        # You can run fastq-dump ERR3958992 to get the fastq file, or for this notebook simply call:
        download_data("https://zenodo.org/records/17207091/files/hifi.ccs.fastq.gz", "data/streptococcus")

    download_SRR_with_progress("SRR24332397", "data/streptococcus", 7.08e9)

