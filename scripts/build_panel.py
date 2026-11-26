"""Build the calibration-panel freeze artefacts (prereg §8 items 4 and 7).

Discharges two freeze-checklist items in docs/prereg_calibration_panel.md:

  §8.7  conserved vs non-conserved position assignment, with its alignment
  §8.4  the n = 48 panel sequence list, as a TSV with a recorded sha256

Three inputs, two of them downloaded and cached under panel/:

  1. vertebrate SS18 orthologues (UniProt)      -> per-position conservation
  2. PDB 7VRB, SMARCA4-SS18 fusion, 2.39 A      -> per-position BRG1 contact
  3. ../eCR_mod_lib/library/module_library.tsv  -> paralogs + off-mechanism controls

Run with --fetch to refresh the cached raw inputs; without it the committed
caches are used, so the freeze is reproducible without network.

Outputs (all under panel/):
  snh_orthologs.tsv     the 70-aa alignment, one row per species
  snh_positions.tsv     per-position conservation x interface x structure call
  calibration_panel_v1.tsv   the panel
"""
import argparse
import csv
import hashlib
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "panel"

# --- fixed parents -----------------------------------------------------------
MOUSE = "MSVAFAAPRQRGKGEITPAAIQKMLDENNHLIQCIMDYQNKGKASECSQYQQILHTNLVYLATIADSNQN"
HUMAN = "MSVAFAAPRQRGKGEITPAAIQKMLDDNNHLIQCIMDSQNKGKTSECSQYQQMLHTNLVYLATIADSNQN"

# --- fixed structural facts, from PDB 7VRB / UniProt Q15532 ------------------
# 7VRB is a single-chain fusion: entity pos 3-44 = SMARCA4 P51532 172-213,
# entity pos 56-143 = SS18 Q15532 14-101. So SS18 residue p is at entity p+42.
SS18_OFFSET = 42
SMARCA4_ENT = range(3, 45)
SS18_ENT = range(56, 144)
HELIX_A1 = range(18, 41)          # UniProt Q15532 HELIX 18-40
HELIX_A2 = range(44, 70)          # UniProt Q15532 HELIX 44-69
INTERFACE_CUTOFF_A = 5.0

# --- pre-registered class thresholds ----------------------------------------
# A position is "conserved" only if it is invariant across every orthologue in
# the alignment, and "variable" only if it falls below VARIABLE_MAX. The gap
# between them is deliberately left unassigned so neither panel class is
# contaminated by borderline positions.
CONSERVED_MIN = 1.00
VARIABLE_MAX = 0.98

SCRAMBLE_SEED = 20261126
N_SCRAMBLES = 2

UNIPROT_QUERY = "gene_exact:SS18 AND taxonomy_id:7742 AND keyword:KW-1185"
UNIPROT_FIELDS = ("accession,reviewed,id,protein_name,organism_name,lineage,"
                  "gene_primary,length,sequence")

# Off-mechanism positive controls (section 4.1, which asks for "a minimal VP16,
# or a top library activator with no BAF link"):
#   ED_VP16_TAD  the prototypical acidic AD. Its library known_interactors are
#                MED25, TBP, TAF4, TAF9 -- Mediator and TFIID, no mSWI/SNF
#                subunit -- so off-mechanism here is recorded, not assumed.
#   ED_VP64      the VP16 tetramer eCR_predictor's own designs use as the opener
#                ED, i.e. the exact domain CLAUDE.md says SNH must not evolve
#                into. Same acidic core as VP16-TAD, so these two are one
#                mechanism, not two.
#   DUX4         an independent strong screen activator (HT-recruit -9.73), not
#                an mSWI/SNF subunit, giving the class a second mechanism.
#
# Selected by module_id, never by name. 1,622 ED names in the library are used
# by more than one record (2,167 records collide), so `name` is not a key:
# "LEUTX_1" alone matches five records, three of them repressor-subtype, and
# keying by name picks ED_TYCKO_2025_RDS_LEUTX_1 -- a repressor. This is defect
# F8a seen from the data side.
OFF_MECHANISM = ["ED_VP16_TAD", "ED_VP64", "ED_TYCKO_2025_ADS_DUX4_1"]
PARALOGS = ["ED_TYCKO_2025_ADS_SS18L1_1", "ED_TYCKO_2025_ADS_SS18L2_1"]

CLADES = ["Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria",
          "Amphibia", "Actinopterygii", "Chondrichthyes", "Cyclostomata",
          "Coelacanthimorpha", "Dipnomorpha"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- fetch -------
def fetch_orthologs(dest: Path) -> None:
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urlencode(
        {"query": UNIPROT_QUERY, "fields": UNIPROT_FIELDS, "format": "tsv"})
    with urlopen(url, timeout=300) as r:
        rel = r.headers.get("X-UniProt-Release", "?")
        dest.write_bytes(r.read())
    print("fetched %s (UniProt release %s)" % (dest.name, rel))


def fetch_cif(dest: Path) -> None:
    with urlopen("https://files.rcsb.org/download/7VRB.cif", timeout=300) as r:
        dest.write_bytes(r.read())
    print("fetched %s" % dest.name)


# ------------------------------------------------------------- alignment -----
def build_alignment(raw: Path) -> list:
    """One ungapped 70-aa SNH per species, aligned to the mouse parent."""
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner(mode="local", open_gap_score=-11,
                                    extend_gap_score=-1)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    rows = list(csv.DictReader(open(raw, encoding="utf-8"), delimiter="\t"))
    by_org = defaultdict(list)
    for r in rows:
        by_org[r["Organism"]].append(r)
    reps = [sorted(v, key=lambda r: (0 if r["Reviewed"] == "reviewed" else 1,
                                     -int(r["Length"]), r["Entry"]))[0]
            for v in by_org.values()]

    kept, drops = [], Counter()
    for r in reps:
        full = r["Sequence"]
        if len(full) < 70:
            drops["entry shorter than 70 aa"] += 1
            continue
        aln = aligner.align(MOUSE, full)[0]
        qa, ta = str(aln[0]), str(aln[1])
        if len(qa.replace("-", "")) < 70:
            # The SNH is N-terminally truncated in the entry, not divergent:
            # these are incomplete gene predictions that begin inside the domain.
            drops["SNH incomplete in entry"] += 1
            continue
        if "-" in qa or "-" in ta:
            drops["gapped alignment"] += 1
            continue
        lin = r["Taxonomic lineage"]
        clade = next((c for c in CLADES if c in lin), "other")
        kept.append({"accession": r["Entry"], "organism": r["Organism"],
                     "reviewed": r["Reviewed"], "clade": clade,
                     "pct_id_to_mouse": round(sum(a == b for a, b in zip(qa, ta)) / 70.0, 4),
                     "snh70": ta})

    print("species considered: %d   aligned: %d" % (len(reps), len(kept)))
    for k, v in drops.most_common():
        print("   dropped, %-28s %d" % (k, v))
    for acc, exp in (("Q62280", MOUSE), ("Q15532", HUMAN)):
        got = next((k["snh70"] for k in kept if k["accession"] == acc), None)
        if got != exp:
            sys.exit("ABORT: %s did not round-trip to its known SNH sequence" % acc)
    print("   parents round-trip exactly: Q62280, Q15532")
    return kept


def conservation(aln: list) -> dict:
    n = len(aln)
    out = {}
    for i in range(1, 71):
        col = Counter(a["snh70"][i - 1] for a in aln)
        mouse_aa = MOUSE[i - 1]
        alt = [(aa, c) for aa, c in col.most_common() if aa != mouse_aa]
        out[i] = {"n": n, "mouse_aa": mouse_aa, "human_aa": HUMAN[i - 1],
                  "frac_identical": col[mouse_aa] / n,
                  "n_distinct": len(col),
                  "top_alt_aa": alt[0][0] if alt else "",
                  "top_alt_frac": (alt[0][1] / n) if alt else 0.0,
                  "counts": ";".join("%s:%d" % kv for kv in col.most_common())}
    return out


# -------------------------------------------------------------- interface ----
def interface(cif: Path) -> dict:
    from Bio.PDB import MMCIFParser

    model = next(iter(MMCIFParser(QUIET=True).get_structure("7vrb", str(cif))))
    chains = {}
    for ch in model:
        d = {r.id[1]: [a for a in r if a.element != "H"]
             for r in ch if r.id[0] == " "}
        if d:
            chains[ch.id] = d

    def mind(aa, bb):
        best = 1e9
        for a in aa:
            ax, ay, az = a.coord
            for b in bb:
                d = (ax - b.coord[0]) ** 2 + (ay - b.coord[1]) ** 2 + (az - b.coord[2]) ** 2
                best = min(best, d)
        return math.sqrt(best)

    intra, inter = defaultdict(lambda: 1e9), defaultdict(lambda: 1e9)
    for cs, ds in chains.items():
        for es, ats in ds.items():
            if es not in SS18_ENT:
                continue
            for cb, db in chains.items():
                for eb, atb in db.items():
                    if eb not in SMARCA4_ENT:
                        continue
                    d = mind(ats, atb)
                    tgt = intra if cs == cb else inter
                    tgt[es] = min(tgt[es], d)

    out = {}
    for ent in sorted(set(intra) | set(inter)):
        p = ent - SS18_OFFSET
        if 1 <= p <= 70:
            out[p] = {"intra": intra.get(ent, float("inf")),
                      "inter": inter.get(ent, float("inf"))}
    n_int = sum(1 for v in out.values() if v["intra"] <= INTERFACE_CUTOFF_A)
    n_ext = sum(1 for v in out.values() if v["inter"] <= INTERFACE_CUTOFF_A)
    print("SNH residues resolved in 7VRB: %d   intra-chain contacts: %d   "
          "inter-chain: %d" % (len(out), n_int, n_ext))
    return out


def element(p: int) -> str:
    if p in HELIX_A1:
        return "helix_a1"
    if p in HELIX_A2:
        return "helix_a2"
    if p in (41, 42, 43):
        return "interhelix_linker"
    return "n_terminal_coil"


def positions(cons: dict, iface: dict) -> list:
    rows = []
    for p in range(1, 71):
        c = cons[p]
        it = iface.get(p)
        resolved = it is not None
        dist = it["intra"] if resolved else None
        is_iface = bool(resolved and dist <= INTERFACE_CUTOFF_A)
        if c["frac_identical"] >= CONSERVED_MIN:
            cons_class = "invariant"
        elif c["frac_identical"] <= VARIABLE_MAX:
            cons_class = "variable"
        else:
            cons_class = "unassigned"
        rows.append({
            "pos": p, "mouse_aa": c["mouse_aa"], "human_aa": c["human_aa"],
            "n_orthologs": c["n"], "frac_identical_to_mouse": round(c["frac_identical"], 4),
            "n_distinct_aa": c["n_distinct"], "top_alt_aa": c["top_alt_aa"],
            "top_alt_frac": round(c["top_alt_frac"], 4),
            "conservation_class": cons_class,
            "resolved_in_7vrb": "yes" if resolved else "no",
            "brg1_min_dist_A": round(dist, 2) if resolved else "",
            "brg1_interface": "yes" if is_iface else ("no" if resolved else "unknown"),
            "structural_element": element(p),
            "aa_counts": c["counts"],
        })
    return rows


# ------------------------------------------------------------------ panel ----
def library_records(lib: Path, module_ids: list) -> dict:
    """Look records up by module_id, which is unique, and prove it is."""
    csv.field_size_limit(10 ** 7)
    want = set(module_ids)
    out, seen = {}, Counter()
    with open(lib, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["type"] != "ED":
                continue
            seen[r["module_id"]] += 1
            if r["module_id"] in want:
                out[r["module_id"]] = r
    ambiguous = [m for m in want if seen[m] > 1]
    if ambiguous:
        sys.exit("ABORT: module_id is not unique for %s" % ", ".join(ambiguous))
    missing = want - set(out)
    if missing:
        sys.exit("ABORT: not in library: %s" % ", ".join(sorted(missing)))
    return out


def mutate(seq: str, pos: int, new: str) -> str:
    return seq[:pos - 1] + new + seq[pos:]


def build_panel(posrows: list, lib: Path) -> list:
    by_pos = {r["pos"]: r for r in posrows}
    panel = []

    def add(pid, cls, in_rho, seq, **kw):
        row = {"panel_id": pid, "class": cls,
               "in_primary_rho": "yes" if in_rho else "no",
               "mutation": kw.get("mutation", ""), "position": kw.get("position", ""),
               "wt_aa": kw.get("wt_aa", ""), "mut_aa": kw.get("mut_aa", ""),
               "sequence_aa": seq, "length_aa": len(seq),
               "frac_identical_to_mouse": kw.get("cons", ""),
               "brg1_min_dist_A": kw.get("dist", ""),
               "brg1_interface": kw.get("iface", ""),
               "structural_element": kw.get("elem", ""),
               "expect_raw_reporter": kw["e_raw"],
               "expect_differential": kw["e_diff"],
               "expect_occupancy": kw["e_occ"],
               "source": kw.get("source", "this repo"),
               "notes": kw.get("notes", "")}
        panel.append(row)

    # 1. parents
    add("EV_SNH_WT_MOUSE", "parent_wt", True, MOUSE,
        e_raw="reference", e_diff="reference", e_occ="reference",
        source="eCR_mod_lib ED_SS18_SNH (Q62280 1-70)",
        notes="anchors both axes; replicate count in prereg section 7.3")
    add("EV_SNH_WT_HUMAN", "parent_ortholog", True, HUMAN,
        e_raw="~parent", e_diff="~parent", e_occ="~parent",
        source="eCR_mod_lib ED_SS18_N70 (Q15532 1-70)",
        notes="natural near-neutral variant; differs at 27/38/44/53, all off-interface")

    # 2. paralogs, trimmed to the SNH boundary
    recs = library_records(lib, PARALOGS + OFF_MECHANISM)
    for mid in PARALOGS:
        r = recs[mid]
        full = r["sequence_aa"]
        add("EV_SNH_%s" % r["name"].replace("_1", ""), "paralog", True, full[:70],
            e_raw="active", e_diff="intermediate", e_occ="intermediate",
            source="eCR_mod_lib %s (%s), trimmed 1-70" % (r["module_id"], r["uniprot_id"]),
            notes="library record is %d aa; trimmed to the SNH boundary for "
                  "panel comparability, so its library HT-recruit score is NOT "
                  "a panel measurement" % len(full))

    # 3. conserved + interface -> alanine scan
    iface_inv = [r["pos"] for r in posrows
                 if r["brg1_interface"] == "yes" and r["conservation_class"] == "invariant"]
    for p in iface_inv:
        r = by_pos[p]
        wt = r["mouse_aa"]
        mut = "G" if wt == "A" else "A"
        add("EV_SNH_%s%d%s" % (wt, p, mut), "conserved_interface", True,
            mutate(MOUSE, p, mut), mutation="%s%d%s" % (wt, p, mut), position=p,
            wt_aa=wt, mut_aa=mut, cons=r["frac_identical_to_mouse"],
            dist=r["brg1_min_dist_A"], iface="yes", elem=r["structural_element"],
            e_raw="reduced", e_diff="reduced to dead", e_occ="reduced to dead",
            notes="glycine used because the parent residue is alanine" if wt == "A" else "")

    # 4. variable + non-interface -> swap to the commonest orthologue alternative
    var_non = [r["pos"] for r in posrows
               if r["conservation_class"] == "variable"
               and r["brg1_interface"] in ("no", "unknown")]
    for p in var_non:
        r = by_pos[p]
        wt, mut = r["mouse_aa"], r["top_alt_aa"]
        note = ("substitutes the commonest alternative residue seen at this "
                "position across the alignment (%.1f%% of orthologues)"
                % (100 * r["top_alt_frac"]))
        if r["brg1_interface"] == "unknown":
            note += "; position unresolved in 7VRB, non-interface inferred from " \
                    "its being N-terminal to both helices"
        add("EV_SNH_%s%d%s" % (wt, p, mut), "nonconserved_noninterface", True,
            mutate(MOUSE, p, mut), mutation="%s%d%s" % (wt, p, mut), position=p,
            wt_aa=wt, mut_aa=mut, cons=r["frac_identical_to_mouse"],
            dist=r["brg1_min_dist_A"], iface=r["brg1_interface"],
            elem=r["structural_element"],
            e_raw="~parent", e_diff="~parent", e_occ="~parent", notes=note)

    # 5. helix-disrupting prolines, mid-helix, invariant but NOT at the interface,
    #    so the effect is loss of the fold rather than deletion of a contact
    helix_pro = [23, 31, 50, 57, 63]
    for p in helix_pro:
        r = by_pos[p]
        if r["conservation_class"] != "invariant" or r["brg1_interface"] != "no":
            sys.exit("ABORT: helix-disrupting position %d is no longer "
                     "invariant/non-interface" % p)
        wt = r["mouse_aa"]
        add("EV_SNH_%s%dP" % (wt, p), "helix_disrupting", True, mutate(MOUSE, p, "P"),
            mutation="%s%dP" % (wt, p), position=p, wt_aa=wt, mut_aa="P",
            cons=r["frac_identical_to_mouse"], dist=r["brg1_min_dist_A"], iface="no",
            elem=r["structural_element"],
            e_raw="dead", e_diff="dead", e_occ="dead",
            notes="mid-helix proline; section 7.4 reads these as the dead control")

    # 6. off-mechanism positive controls
    for mid in OFF_MECHANISM:
        r = recs[mid]
        slug = "".join(c if c.isalnum() else "_"
                       for c in r["name"].replace("_1", "")).strip("_")
        bits = ["strong library activator, not an mSWI/SNF subunit"]
        if r["known_interactors"]:
            bits.append("library known_interactors: %s" % r["known_interactors"])
        if r["quantitative_metric"]:
            bits.append("%s = %.3f" % (r["quantitative_metric_label"],
                                       float(r["quantitative_metric"])))
        bits.append("off-mechanism is the premise section 7.2 tests; if these do "
                    "not separate, the panel is void by that control")
        add("EV_OFFMECH_%s" % slug, "off_mechanism", True, r["sequence_aa"],
            e_raw="high", e_diff="~zero", e_occ="background",
            source="eCR_mod_lib %s%s" % (r["module_id"],
                                         " (%s)" % r["uniprot_id"]
                                         if r["uniprot_id"] else ""),
            notes="; ".join(bits))

    # 7. floor, excluded from the primary rho per section 5
    add("EV_FLOOR_NO_EFFECTOR", "floor", False, "",
        e_raw="floor", e_diff="floor", e_occ="background",
        notes="tethering domain only, no effector")
    rng = random.Random(SCRAMBLE_SEED)
    for i in range(1, N_SCRAMBLES + 1):
        chars = list(MOUSE)
        rng.shuffle(chars)
        add("EV_FLOOR_SCRAMBLE_%d" % i, "floor", False, "".join(chars),
            e_raw="floor", e_diff="floor", e_occ="background",
            notes="composition-matched shuffle of the mouse parent, "
                  "random.Random(%d) draw %d" % (SCRAMBLE_SEED, i))
    return panel


def check_names(panel: list, lib: Path) -> None:
    """CLAUDE.md inherited constraint: a name that collides with the library is
    silently dropped by eCR_predictor's fuse.py _load_eds (defect F8a)."""
    csv.field_size_limit(10 ** 7)
    taken = set()
    with open(lib, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            taken.add(r["module_id"])
            taken.add(r["name"])
    ids = [r["panel_id"] for r in panel]
    dupes = [k for k, v in Counter(ids).items() if v > 1]
    clash = sorted(set(ids) & taken)
    if dupes:
        sys.exit("ABORT: duplicate panel_id: %s" % ", ".join(dupes))
    if clash:
        sys.exit("ABORT: panel_id collides with eCR_mod_lib: %s" % ", ".join(clash))
    print("panel ids: %d, all unique and none colliding with the library" % len(ids))


def write_tsv(path: Path, rows: list) -> None:
    # LF, explicitly: these files' sha256 values are recorded in the
    # pre-registration, so their bytes must not depend on the platform that
    # wrote them. .gitattributes stops git rewriting them on checkout.
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("wrote %-34s %3d rows  sha256 %s" % (path.name, len(rows), sha256(path)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true",
                    help="refresh the cached UniProt and RCSB inputs")
    ap.add_argument("--library", default=str(ROOT.parent / "eCR_mod_lib" /
                                            "library" / "module_library.tsv"))
    args = ap.parse_args()

    raw = PANEL / "ss18_orthologs_raw.tsv"
    cif = PANEL / "7vrb.cif"
    if args.fetch or not raw.exists():
        fetch_orthologs(raw)
    if args.fetch or not cif.exists():
        fetch_cif(cif)
    lib = Path(args.library)
    if not lib.exists():
        sys.exit("ABORT: library TSV not found at %s" % lib)

    print("\n-- inputs")
    for p in (raw, cif, lib):
        print("   %-30s sha256 %s" % (p.name, sha256(p)))

    print("\n-- alignment")
    aln = build_alignment(raw)
    clades = Counter(a["clade"] for a in aln)
    print("   clades: " + ", ".join("%s %d" % kv for kv in clades.most_common()))
    write_tsv(PANEL / "snh_orthologs.tsv",
              sorted(aln, key=lambda a: (-a["pct_id_to_mouse"], a["organism"])))

    print("\n-- interface")
    iface = interface(cif)

    print("\n-- positions")
    posrows = positions(conservation(aln), iface)
    cc = Counter(r["conservation_class"] for r in posrows)
    print("   conservation: " + ", ".join("%s %d" % kv for kv in cc.most_common()))
    print("   interface: " + ", ".join(
        "%s %d" % kv for kv in Counter(r["brg1_interface"] for r in posrows).most_common()))
    inv_if = [r["pos"] for r in posrows if r["brg1_interface"] == "yes"]
    print("   every interface residue invariant: %s"
          % all(by["conservation_class"] == "invariant"
                for by in posrows if by["pos"] in inv_if))
    write_tsv(PANEL / "snh_positions.tsv", posrows)

    print("\n-- panel")
    panel = build_panel(posrows, lib)
    check_names(panel, lib)
    byc = Counter(r["class"] for r in panel)
    for c, n in byc.most_common():
        print("   %-28s %2d" % (c, n))
    n_rho = sum(1 for r in panel if r["in_primary_rho"] == "yes")
    print("   total constructs %d, of which %d enter the primary rho" % (len(panel), n_rho))
    if n_rho < 48:
        sys.exit("ABORT: prereg section 6 sizes the decision rule at n >= 48, got %d"
                 % n_rho)
    print("   n = %d, at or above the section 6 power floor of 48" % n_rho)
    write_tsv(PANEL / "calibration_panel_v1.tsv", panel)


if __name__ == "__main__":
    main()
