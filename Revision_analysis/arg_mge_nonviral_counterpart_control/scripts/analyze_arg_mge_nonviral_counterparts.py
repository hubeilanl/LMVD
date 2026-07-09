#!/usr/bin/env python3
"""Build ARG-MGE non-viral counterpart evidence tables."""
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

TARGET_PAIRS_RAW = [
    ("tnp", "lnuC"),
    ("tndX", "AAC(6')-Im"),
    ("tndX", "APH(2'')-IIa"),
    ("tnpR_2", "ACI-1"),
]


def norm(name):
    if pd.isna(name):
        return ""
    value = str(name).strip()
    value = value.replace("’", "'").replace("‘", "'")
    value = value.replace("“", '"').replace("”", '"')
    value = value.replace('"', "'")
    return re.sub(r"\s+", " ", value)


TARGET_PAIRS = [(norm(m), norm(a)) for m, a in TARGET_PAIRS_RAW]


def bounds(start, stop):
    return min(start, stop), max(start, stop)


def interval_distance(a1, a2, b1, b2):
    a_start, a_end = bounds(a1, a2)
    b_start, b_end = bounds(b1, b2)
    if a_start <= b_end and b_start <= a_end:
        return 0
    return b_start - a_end if a_end < b_start else a_start - b_end


def is_similar(vdist, ndist):
    diff = abs(ndist - vdist)
    if diff <= 1000:
        return True
    return (diff / max(vdist, 1)) <= 0.20


s1 = pd.read_csv(DATA / "target_pair_viral_arg_annotations.tsv", sep="\t")
s2 = pd.read_csv(DATA / "target_pair_viral_mge_annotations.tsv", sep="\t")
s3 = pd.read_csv(DATA / "target_pair_all_contig_arg_annotations.tsv", sep="\t")
s4 = pd.read_csv(DATA / "target_pair_all_contig_mge_annotations.tsv", sep="\t")

s1.columns = [
    "sample", "original_contig", "viral_contig", "arg_start", "arg_stop",
    "arg_orient", "bitscore", "arg_name", "aro", "drug_class",
    "mechanism", "family",
]
s2.columns = [
    "sample", "original_contig", "viral_contig", "mobileog",
    "mge_start", "mge_stop", "mge_name",
]
s3.columns = ["sample", "contig", "arg_start", "arg_stop", "arg_orient", "arg_name", "aro"]
s4.columns = ["sample", "contig", "mobileog", "mge_start", "mge_stop", "mge_name"]

s1["arg_norm"] = s1["arg_name"].map(norm)
s2["mge_norm"] = s2["mge_name"].map(norm)
s3["arg_norm"] = s3["arg_name"].map(norm)
s4["mge_norm"] = s4["mge_name"].map(norm)

viral = pd.merge(
    s1[["sample", "original_contig", "viral_contig", "arg_start", "arg_stop", "arg_orient", "arg_name", "arg_norm"]],
    s2[["sample", "original_contig", "viral_contig", "mge_start", "mge_stop", "mge_name", "mge_norm"]],
    on=["sample", "original_contig", "viral_contig"],
    how="inner",
)

mask = pd.Series(False, index=viral.index)
for mge_name, arg_name in TARGET_PAIRS:
    mask |= (viral["mge_norm"] == mge_name) & (viral["arg_norm"] == arg_name)
viral = viral[mask].copy()
viral["pair_name"] = viral["mge_norm"] + " + " + viral["arg_norm"]
viral["distance"] = viral.apply(
    lambda r: interval_distance(r.arg_start, r.arg_stop, r.mge_start, r.mge_stop),
    axis=1,
)
viral["within_5kb"] = viral["distance"].le(5000).map({True: "yes", False: "no"})
viral["within_10kb"] = viral["distance"].le(10000).map({True: "yes", False: "no"})
viral = viral[viral["distance"] <= 10000].reset_index(drop=True)

allc = pd.merge(
    s3[["sample", "contig", "arg_start", "arg_stop", "arg_orient", "arg_name", "arg_norm"]],
    s4[["sample", "contig", "mge_start", "mge_stop", "mge_name", "mge_norm"]],
    on=["sample", "contig"],
    how="inner",
)
mask = pd.Series(False, index=allc.index)
for mge_name, arg_name in TARGET_PAIRS:
    mask |= (allc["mge_norm"] == mge_name) & (allc["arg_norm"] == arg_name)
allc = allc[mask].copy()
allc["pair_name"] = allc["mge_norm"] + " + " + allc["arg_norm"]
allc["distance"] = allc.apply(
    lambda r: interval_distance(r.arg_start, r.arg_stop, r.mge_start, r.mge_stop),
    axis=1,
)
allc["within_5kb"] = allc["distance"].le(5000).map({True: "yes", False: "no"})
allc["within_10kb"] = allc["distance"].le(10000).map({True: "yes", False: "no"})

viral_original_contigs = {}
for _, row in pd.concat([s1[["sample", "original_contig"]], s2[["sample", "original_contig"]]]).iterrows():
    viral_original_contigs.setdefault(row["sample"], set()).add(row["original_contig"])


def is_viral_contig(row):
    return row["contig"] in viral_original_contigs.get(row["sample"], set())


nonviral = allc[~allc.apply(is_viral_contig, axis=1)].copy()

main_rows = []
for i, vu in viral.iterrows():
    nv_same = nonviral[
        (nonviral["sample"] == vu["sample"]) &
        (nonviral["pair_name"] == vu["pair_name"])
    ].copy()
    cp10 = nv_same[nv_same["distance"] <= 10000]
    cp5 = nv_same[nv_same["distance"] <= 5000]
    cp_sim = nv_same[nv_same["distance"].map(lambda d: is_similar(vu["distance"], d))]
    closest = nv_same.sort_values("distance").head(1)
    closest_row = closest.iloc[0] if len(closest) else None
    main_rows.append({
        "unit_id": f"unit_{i + 1:03d}",
        "sample": vu["sample"],
        "original_contig": vu["original_contig"],
        "viral_contig": vu["viral_contig"],
        "pair_name": vu["pair_name"],
        "ARG_name": vu["arg_name"],
        "MGE_name": vu["mge_name"],
        "ARG_start_viral": vu["arg_start"],
        "ARG_stop_viral": vu["arg_stop"],
        "ARG_orientation_viral": vu["arg_orient"],
        "MGE_start_viral": vu["mge_start"],
        "MGE_stop_viral": vu["mge_stop"],
        "ARG_MGE_distance_viral_bp": vu["distance"],
        "within_5kb_viral": vu["within_5kb"],
        "within_10kb_viral": vu["within_10kb"],
        "same_sample_nonviral_counterpart_within_10kb": "yes" if len(cp10) else "no",
        "same_sample_nonviral_counterpart_within_5kb": "yes" if len(cp5) else "no",
        "same_sample_nonviral_counterpart_similar_distance": "yes" if len(cp_sim) else "no",
        "number_of_same_sample_nonviral_counterparts_within_10kb": len(cp10),
        "number_of_same_sample_nonviral_counterparts_within_5kb": len(cp5),
        "number_of_same_sample_nonviral_counterparts_similar_distance": len(cp_sim),
        "closest_same_sample_nonviral_contig": closest_row["contig"] if closest_row is not None else "NA",
        "closest_same_sample_nonviral_distance_bp": closest_row["distance"] if closest_row is not None else "NA",
        "control_interpretation": (
            "Same-sample non-viral counterpart detected"
            if len(cp10)
            else "No same-sample non-viral counterpart detected for the same ARG-MGE pair within 10 kb"
        ),
    })

main = pd.DataFrame(main_rows)
main.to_csv(RESULTS / "viral_arg_mge_units_with_nonviral_counterpart_status.tsv",
            sep="\t", index=False)

nv_rows = []
for _, row in nonviral.iterrows():
    viral_same = viral[
        (viral["sample"] == row["sample"]) &
        (viral["pair_name"] == row["pair_name"])
    ]
    nv_rows.append({
        "sample": row["sample"],
        "nonviral_contig": row["contig"],
        "pair_name": row["pair_name"],
        "ARG_name": row["arg_name"],
        "MGE_name": row["mge_name"],
        "ARG_MGE_distance_nonviral_bp": row["distance"],
        "same_sample_counterpart_to_viral_unit": "yes" if len(viral_same) else "no",
        "number_of_same_sample_viral_units": len(viral_same),
    })
pd.DataFrame(nv_rows).to_csv(
    RESULTS / "same_sample_nonviral_candidate_pairs.tsv",
    sep="\t", index=False,
)

by_pair = []
for pair_name, sub in main.groupby("pair_name"):
    by_pair.append({
        "pair_name": pair_name,
        "viral_unit_count": len(sub),
        "same_sample_counterpart_10kb_count": (sub["same_sample_nonviral_counterpart_within_10kb"] == "yes").sum(),
        "same_sample_counterpart_10kb_percent": (sub["same_sample_nonviral_counterpart_within_10kb"] == "yes").mean() * 100,
    })
by_pair = pd.DataFrame(by_pair)
by_pair.to_csv(RESULTS / "nonviral_counterpart_summary_by_pair.tsv",
               sep="\t", index=False)

overall = pd.DataFrame([{
    "viral_unit_count": len(main),
    "sample_count": main["sample"].nunique(),
    "pair_count": main["pair_name"].nunique(),
    "same_sample_counterpart_10kb_count": (main["same_sample_nonviral_counterpart_within_10kb"] == "yes").sum(),
    "same_sample_counterpart_10kb_percent": (main["same_sample_nonviral_counterpart_within_10kb"] == "yes").mean() * 100,
    "nonviral_candidate_pair_count": len(nonviral),
}])
overall.to_csv(RESULTS / "nonviral_counterpart_summary_overall.tsv",
               sep="\t", index=False)

summary = [
    "ARG-MGE same-sample non-viral counterpart control",
    f"Viral target ARG-MGE units evaluated: {len(main)}",
    f"Samples represented: {main['sample'].nunique()}",
    f"Same-sample non-viral counterparts within 10 kb: {(main['same_sample_nonviral_counterpart_within_10kb'] == 'yes').sum()}",
    f"Non-viral candidate pairs for target combinations: {len(nonviral)}",
]
(RESULTS / "arg_mge_nonviral_counterpart_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))
