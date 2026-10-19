"""_arm_vocab.py — the HH21 module positive control's arm names and their test directions, in ONE place.

`arm` itself is minted by the promoter scan (`grn_memory_node_scan.py`, beside the flags that
define it) and travels in `response_peakgroups.tsv` / `hh21_hh1716_response_counts.tsv`. What
lives here is the small vocabulary the two Python consumers — the module test (which scores the
arms) and the panel (which plots them and writes the label-swap stats) — both need to agree on:
which arms are SCORED, which is held out, and each arm's alternative hypothesis.

They used to declare it twice. Name drift failed loud (a missing module raises), but the
DIRECTION map did not: flipping PITX2ab_ctrl to one-sided in one file would have produced
`hh21_hh1716_module_stats.tsv` and `hh21_hh1716_module_labelswap_stats.tsv` disagreeing about
`alt` for the same arm, with no error — the two-p's-on-one-count class of inconsistency that
 removed elsewhere.

The R model (`hh21_hh1716_nb_glmm.R`) re-declares the names on its side of the language
boundary and validates them against the data, which is the right shape for a cross-language
contract: it stops on any unexpected arm value.
"""

L_ARM = "L-up"                # response arm: promoters of genes Left-up in the HH21 RNA
R_ARM = "R-up"                # response arm: promoters of genes Right-up in the HH21 RNA
REF_ARM = "PITX2c_ref"        # the PITX2c promoter PEAK (fig03's headline element)
AB_ARM = "PITX2ab_ctrl"       # PITX2a/b promoters — bilateral within-locus control
LOCUS_ARM = "PITX2_locus_other"   # PITX2's remaining cis-elements, HELD OUT of the response arms

# Scored arms, in panel order. LOCUS_ARM is deliberately absent: it is carried in the counts
# table for provenance but never scored, so the module tests genes other than PITX2.
MODULES = [L_ARM, R_ARM, REF_ARM, AB_ARM]

# One-sided where the biology predicts a direction; two-sided for the control, which predicts
# none — so its star is never a one-sided pass.
ALT = {
    L_ARM: "greater",         # L-up genes' promoters should be Left-open
    R_ARM: "less",            # R-up genes' promoters should be Right-open
    REF_ARM: "greater",       # PITX2c is the Left-biased reference
    AB_ARM: "two-sided",      # a/b are bilateral; no direction is predicted
}
ALT_DIR = "greater"           # L-up-minus-R-up directionality: Left-ward is the prediction

assert set(ALT) == set(MODULES), "every scored arm needs an alternative"
