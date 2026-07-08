"""
Pull Allen Mouse Brain Atlas ISH expression for genes underlying the three
circuits at the center of the Ganguly lab's work:

  1. Sleep spindle / slow-oscillation generators
       Cacna1g (Cav3.1 — thalamic T-type Ca2+, drives spindles)
       Cacna1h (Cav3.2)
       Hcn1    (Ih current — pacemaker for thalamocortical rhythm)
       Kcnq2   (M-current — cortical excitability)

  2. GABAergic consolidation (sleep + stroke recovery)
       Gad1    (GAD67 — primary GABA synthesis)
       Gad2    (GAD65)
       Pvalb   (parvalbumin interneurons — fast-spiking, spindle-coupled)
       Sst     (somatostatin interneurons)

  3. Corticostriatal motor learning
       Drd1    (D1 receptor — direct pathway, reward/consolidation)
       Drd2    (D2 receptor — indirect pathway)
       Foxp1   (MSN marker — striatal output neurons)
       Slc17a7 (VGLUT1 — corticostriatal excitatory input)

Structures: primary motor cortex (MOp 985), anterior cingulate/PFC (ACA 31),
            striatum (STR 672), thalamus (TH 549), hippocampus (HPF 1089),
            substantia nigra (SN 374), brainstem (BS 313)
"""
import requests, csv, time
from pathlib import Path

OUT  = Path(__file__).parent
BASE = "http://api.brain-map.org/api/v2/data/query.json"

GENES = {
    "spindle_osc": ["Cacna1g", "Cacna1h", "Hcn1", "Kcnq2"],
    "gaba":        ["Gad1", "Gad2", "Pvalb", "Sst"],
    "corticostr":  ["Drd1", "Drd2", "Foxp1", "Slc17a7"],
}
ALL_GENES = [g for grp in GENES.values() for g in grp]

STRUCTURES = {
    "Motor_ctx":   985,
    "PFC":         31,
    "Striatum":    672,
    "Thalamus":    549,
    "Hippocampus": 1089,
    "Subst_nigra": 374,
    "Midbrain":    313,
}

def get_dataset(gene):
    url = (f"{BASE}?criteria=model::SectionDataSet,"
           f"rma::criteria,genes[acronym$eq'{gene}'],"
           f"[failed$eqfalse],"
           f"products[abbreviation$eqMouse]&num_rows=5")
    r = requests.get(url, timeout=20).json()
    if not r.get("msg"):
        return None
    for row in r["msg"]:
        if row.get("plane_of_section_id") == 2:
            return row["id"]
    return r["msg"][0]["id"] if r["msg"] else None

def get_expression(dataset_id, struct_ids):
    id_str = ",".join(str(i) for i in struct_ids)
    url = (f"{BASE}?criteria=model::StructureUnionize,"
           f"rma::criteria,[section_data_set_id$eq{dataset_id}],"
           f"[structure_id$in{id_str}]&num_rows=50")
    r = requests.get(url, timeout=20).json()
    return {row["structure_id"]: row.get("expression_energy", 0.0) or 0.0
            for row in r.get("msg", [])}

struct_ids = list(STRUCTURES.values())
id_to_name = {v: k for k, v in STRUCTURES.items()}
gene_to_grp = {g: grp for grp, gs in GENES.items() for g in gs}

rows = []
for gene in ALL_GENES:
    print(f"  {gene}...", end=" ", flush=True)
    ds = get_dataset(gene)
    if ds is None:
        print("no dataset"); continue
    expr = get_expression(ds, struct_ids)
    row = {"gene": gene, "group": gene_to_grp[gene], "dataset_id": ds}
    for sid, name in id_to_name.items():
        row[name] = round(expr.get(sid, 0.0), 4)
    rows.append(row)
    print(f"ds={ds} ✓")
    time.sleep(0.3)

fieldnames = ["gene","group","dataset_id"] + list(STRUCTURES.keys())
with open(OUT / "expression.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"\nSaved expression.tsv ({len(rows)} genes)")
for r in rows:
    vals = {k: r[k] for k in STRUCTURES}
    top  = max(vals, key=vals.get)
    print(f"  {r['gene']:12s}  peak={top} ({vals[top]:.3f})")
