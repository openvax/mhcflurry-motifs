"""
Generate allele pages and the main page

"""
import sys
import argparse
import os
import shutil
import jinja2

import pandas
import tqdm
parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "--motif-artifacts",
    metavar="DIR",
    default="current/analysis_predictor_info/motifs",
    help="Motifs artifacts DIR. Default: %(default)s.")
parser.add_argument(
    "--templates",
    metavar="DIR",
    default="templates/",
    help="Templates dir. Default: %(default)s.")
parser.add_argument(
    "--version-file",
    metavar="TXT",
    default="current/MHCFLURRY_VERSION.txt",
    help="Version string. Default: %(default)s.")
parser.add_argument(
    "--max-alleles",
    metavar="N",
    type=int,
    help="Use only N alleles. For debugging.")
parser.add_argument(
    "--out",
    default="docs/",
    metavar="DIR",
    help="Out DIR. Default: %(default)s.")


def page_name(allele):
    name = "%s.html" % (
        allele.replace("*", "-").replace(":", "-"))
    return name


def run():
    args = parser.parse_args(sys.argv[1:])

    if not os.path.exists(args.out):
        os.mkdir(args.out)
        print("Created: ", args.out)

    motif_artifacts_df = pandas.read_csv(
        os.path.join(args.motif_artifacts, "artifacts.csv"))
    motif_artifacts_df["page"] = motif_artifacts_df.allele.map(page_name)
    motif_artifacts_df["sort_key"] = ~motif_artifacts_df.allele.str.startswith("HLA")
    motif_artifacts_df = motif_artifacts_df.sort_values(["sort_key", "allele"])

    version = open(args.version_file, "r").read()

    if args.max_alleles:
        motif_artifacts_df = motif_artifacts_df.head(args.max_alleles)

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(args.templates))
    index_template = env.get_template('index.html')

    out_filename = os.path.join(args.out, 'index.html')
    with open(out_filename, 'w') as fh:
        fh.write(
            index_template.render(
                page_name="home",
                version=version,
                motif_artifacts_df=motif_artifacts_df))

    print("Wrote %s" % out_filename)

    print("Writing allele pages.")
    allele_template = env.get_template('allele.html')
    for _, row in tqdm.tqdm(
            motif_artifacts_df.iterrows(), total=len(motif_artifacts_df)):
        out_filename = os.path.join(args.out, row.page)
        d = row.to_dict()
        with open(out_filename, 'w') as fh:
            fh.write(
                allele_template.render(
                    page_name="allele",
                    version=version,
                    motif_artifacts_df=motif_artifacts_df,
                    **d))

        for item in ['logo_filename', 'length_distribution_filename']:
            # Make a hard link
            if os.path.exists(os.path.join(args.out, row[item])):
                os.unlink(os.path.join(args.out, row[item]))
            os.link(
                os.path.join(args.motif_artifacts, "artifacts", row[item]),
                os.path.join(args.out, row[item]))

    print("Done.")

if __name__ == '__main__':
    run()
