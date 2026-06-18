"""
DZ-twinning polygenic predisposition score — web app.

Loads the model bundle produced by the analysis notebook and computes a
literature-weighted polygenic score (PRS) from four maternal DZ-twinning loci.

IMPORTANT (and surfaced in the UI): this is a RESEARCH predisposition score,
not a clinical or reproductive predictor. The leakage-corrected analysis showed
these loci add no predictive signal beyond genome-wide ancestry.
"""
import math
import os
import joblib
from flask import Flask, render_template, request

BASE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = joblib.load(os.path.join(BASE, "twinning_model_bundle.joblib"))

SCORE_TABLE = BUNDLE["score_table"]          # rsid, gene, effect_allele, OR, beta
BETAS = BUNDLE["betas"]                       # {rsid: beta}
SNPS = BUNDLE["snps"]                         # ordered list of rsids
PRS_MEAN = float(BUNDLE["prs_mean"])          # reference mean (1000 Genomes)
PRS_STD = float(BUNDLE["prs_std"])            # reference SD (1000 Genomes)
DATA_SOURCE = BUNDLE.get("data_source", "1000 Genomes")

# Per-locus metadata for the form (effect allele = twinning-increasing allele)
LOCI = []
for rsid in SNPS:
    row = SCORE_TABLE[SCORE_TABLE["rsid"] == rsid].iloc[0]
    LOCI.append({
        "rsid": rsid,
        "gene": str(row["gene"]),
        "effect_allele": str(row["effect_allele"]),
        "or": float(row["OR"]),
        "beta": float(BETAS[rsid]),
    })

app = Flask(__name__)


def normal_cdf(z):
    """Standard-normal CDF via erf (no scipy dependency)."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def score_genotypes(dosages):
    """dosages: {rsid: 0|1|2 copies of the effect allele} -> result dict."""
    contributions = []
    prs = 0.0
    for locus in LOCI:
        d = int(dosages[locus["rsid"]])
        contrib = d * locus["beta"]
        prs += contrib
        contributions.append({**locus, "dosage": d, "contribution": contrib})
    z = (prs - PRS_MEAN) / PRS_STD if PRS_STD else 0.0
    percentile = normal_cdf(z) * 100.0
    return {
        "prs": prs,
        "z": z,
        "percentile": percentile,
        "contributions": contributions,
    }


@app.route("/", methods=["GET"])
def index():
    return render_template("index.html", loci=LOCI, result=None,
                           data_source=DATA_SOURCE)


@app.route("/score", methods=["POST"])
def score():
    dosages, errors = {}, []
    for locus in LOCI:
        raw = request.form.get(locus["rsid"], "")
        if raw not in {"0", "1", "2"}:
            errors.append(f"Choose 0, 1, or 2 copies for {locus['gene']} "
                          f"({locus['rsid']}).")
        else:
            dosages[locus["rsid"]] = raw
    if errors:
        return render_template("index.html", loci=LOCI, result=None,
                               errors=errors, data_source=DATA_SOURCE), 400
    result = score_genotypes(dosages)
    return render_template("index.html", loci=LOCI, result=result,
                           submitted=dosages, data_source=DATA_SOURCE)


if __name__ == "__main__":
    # http://127.0.0.1:5000
    app.run(debug=True)
