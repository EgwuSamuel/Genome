# DZ-Twinning Predisposition Score — Flask app

A small web app that computes a literature-weighted polygenic score (PRS) for
maternal dizygotic (DZ) twinning from four loci (FSHB, SMAD3, GNRH1, ZFPM1) and
places it against the 1000 Genomes reference distribution.

**This is a research predisposition score, not a clinical or reproductive
predictor.** The analysis behind it showed these loci add no predictive signal
beyond genome-wide ancestry; the UI states this prominently.

## Files
- `app.py` — Flask backend (loads the model bundle, computes PRS, percentile, per-locus contributions)
- `templates/index.html` — single-page UI (form + result + disclaimer)
- `twinning_model_bundle.joblib` — the trained bundle exported from the analysis notebook
- `requirements.txt` — dependencies

## Run locally
```bash
cd twinning_flask_app
python -m venv venv && source venv/bin/activate      # Windows: venv\Scripts\activate
pip install -r requirements.txt
python app.py
```
Then open http://127.0.0.1:5000 in your browser.

## How to use
For each locus, choose how many copies (0, 1, or 2) of the **effect allele**
(the twinning-increasing allele, shown next to each gene) the person carries.
Click **Compute score** to see the PRS, its z-score, an approximate percentile
versus the reference panel, and each locus's contribution.

## Updating the model
Re-run the analysis notebook to regenerate `twinning_model_bundle.joblib`
(same filename) and drop it in this folder. The app reads `score_table`,
`betas`, `snps`, `prs_mean`, and `prs_std` from the bundle, so adding loci in
the notebook automatically updates the form here.

## Deploying (optional)
Any host that runs Python works (Render, Railway, PythonAnywhere, a VPS with
gunicorn). For production, run behind gunicorn, e.g.:
```bash
pip install gunicorn
gunicorn app:app --bind 0.0.0.0:8000
```
