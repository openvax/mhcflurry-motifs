"""
Generate allele pages and the main page

"""
import sys
import argparse
import os
import shutil
import glob
import markdown
import jinja2

import pandas
import tqdm
parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "--analysis-predictor-info",
    metavar="DIR",
    default="current/analysis_predictor_info",
    help="Analysis predictor info download. Default: %(default)s.")
parser.add_argument(
    "--models-class1-pan",
    metavar="DIR",
    default="current/models_class1_pan",
    help="BA models download. Default: %(default)s.")
parser.add_argument(
    "--templates",
    metavar="DIR",
    default="templates/",
    help="Templates dir. Default: %(default)s.")
parser.add_argument(
    "--markdown",
    metavar="DIR",
    default="markdown/",
    help="Markdown dir. Default: %(default)s.")
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


PAGES = [
    'index.html',
    'about.html',
    'datasets.html',
    'contact.html',
]


ADDITIONAL_MODELS_CLASS1_PAN_FILES = {
    'models.combined/train_data.csv.bz2': 'mhcflurry.ba.train_data.csv.bz2',
    'models.combined/length_distributions.csv.bz2': 'mhcflurry.ba.length_distributions.csv.bz2',
    'models.combined/allele_sequences.csv': 'mhcflurry.allele_sequences.csv',
    'models.combined/frequency_matrices.csv.bz2': 'mhcflurry.ba.frequency_matrices.csv.bz2',
}


def run():
    args = parser.parse_args(sys.argv[1:])

    if not os.path.exists(args.out):
        os.mkdir(args.out)
        print("Created: ", args.out)

    motif_artifacts_df = pandas.read_csv(
        os.path.join(args.analysis_predictor_info, "motifs", "artifacts.csv"))
    motif_artifacts_df["page"] = motif_artifacts_df.allele.map(page_name)

    # Expand motif_artifacts_df to include a row for every allele, even those
    # that have the same allele sequence as another.
    allele_sequences = pandas.read_csv(
        os.path.join(
            args.models_class1_pan, 'models.combined', 'allele_sequences.csv'),
        index_col=0)
    allele_sequences['allele2'] = allele_sequences.index
    motif_artifacts_df["sequence"] = motif_artifacts_df.allele.map(
        allele_sequences.sequence)
    motif_artifacts_df = pandas.merge(
        motif_artifacts_df,
        allele_sequences,
        on="sequence",
        how="outer")
    motif_artifacts_df["base_allele"] = motif_artifacts_df["allele"]
    motif_artifacts_df["allele"] = motif_artifacts_df["allele2"]
    motif_artifacts_df["redundant"] = (
            motif_artifacts_df["base_allele"] != motif_artifacts_df["allele"])

    equivalent_alleles = {}
    sequence_to_alleles = allele_sequences.groupby("sequence").allele2.unique()
    for equivalence_set in sequence_to_alleles.values:
        for allele in equivalence_set:
            equivalent_alleles[allele] = sorted([
                a for a in equivalence_set if a != allele
            ])

    motif_artifacts_df["sort_key"] = (
            (~motif_artifacts_df.allele.str.startswith("HLA")))
    motif_artifacts_df = motif_artifacts_df.sort_values(
        ["sort_key", "redundant", "allele"])

    if args.max_alleles:
        motif_artifacts_df = motif_artifacts_df.head(args.max_alleles)

    rendered_markdown = {}
    if args.markdown:
        for f in glob.glob(os.path.join(args.markdown, "*.md")):
            with open(f, "r", encoding="utf-8") as fd:
                text = fd.read()
            html = markdown.markdown(text)
            name = os.path.basename(f).replace(".md", "")
            rendered_markdown[name] = html

    common_variables = {
        'version': open(args.version_file, "r").read(),
        'motif_artifacts_df': motif_artifacts_df,
        'equivalent_alleles': equivalent_alleles,
        'markdown': rendered_markdown,
    }

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(args.templates))

    print("Writing main pages.")
    for page in PAGES:
        page_template = env.get_template(page)
        out_filename = os.path.join(args.out, page)
        with open(out_filename, 'w') as fh:
            fh.write(
                page_template.render(
                    page_name=page.replace(".html", ""),
                    **common_variables))
        print("Wrote %s" % out_filename)

    print("Writing allele pages.")
    allele_template = env.get_template('allele.html')
    non_redundant = motif_artifacts_df.loc[~motif_artifacts_df.redundant]
    for _, row in tqdm.tqdm(
            non_redundant.iterrows(),
            total=len(non_redundant)):
        out_filename = os.path.join(args.out, row.page)
        d = row.to_dict()
        d.update(common_variables)
        with open(out_filename, 'w') as fh:
            fh.write(
                allele_template.render(
                    page_name="allele",
                    **d))

        for item in ['logo_filename', 'length_distribution_filename']:
            # Make a hard link
            if os.path.exists(os.path.join(args.out, row[item])):
                os.unlink(os.path.join(args.out, row[item]))
            os.link(
                os.path.join(args.analysis_predictor_info, "motifs", "artifacts", row[item]),
                os.path.join(args.out, row[item]))

    # Additional hard links
    for (source, dest) in ADDITIONAL_MODELS_CLASS1_PAN_FILES.items():
        # Make a hard link
        if os.path.exists(os.path.join(args.out, dest)):
            os.unlink(os.path.join(args.out, dest))
        os.link(
            os.path.join(args.models_class1_pan, source),
            os.path.join(args.out, dest))
        print("Linked", os.path.join(args.out, dest))

    print("Done.")


if __name__ == '__main__':
    run()
