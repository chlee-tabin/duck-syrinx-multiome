# Convenience targets for the duck-syrinx-multiome repo.
#
# The primary developer workflow is through the cluster submission wrapper
# and the figures/ directory; this Makefile exposes a handful of
# repeatable checks.

.PHONY: check check-derived-data check-portable-paths check-cluster-python check-lr-binom check-module-inputs check-same-side-duplicates check-locus-obs install-hooks lint-scripts

# ── Derived-data freshness check ───────────────────────────────────────
# Runs the local validator — same one that pre-commit runs.
check-derived-data:
	python3 scripts/verify_derived_data.py

# ── Portable-path guardrail ────────────────────────────────────────────
# Fails on a hardcoded single-user root (/n/scratch/users/<u>/..., /n/standby/...,
# /home/<u>) outside the sanctioned portable idioms.  pre-commit runs the same
# validator on staged files only; this target sweeps the whole repo.
# See docs/audit_hardcoded_paths.md.
check-portable-paths:
	python3 scripts/verify_portable_paths.py --all

# Cluster-executed Python must compile on the CLUSTER's interpreter (3.11), not just this
# workstation's. A 3.12-only f-string once passed local py_compile, synced, queued and died.
check-cluster-python:
	python3 scripts/verify_cluster_python.py

# ── Estimator self-check ───────────────────────────────────────────────
# pipeline/_lr_binom.py is the one depth-adjusted two-sided binomial behind g2c, h4
# and h8; its __main__ asserts direction balance under a true null plus the guard
# branches. Nothing else exercised it (~7 s).
# numpy/scipy live in the analysis env, not the system python the other targets use;
# override with `make check-lr-binom PYTHON_SCI="python3"` if yours has them.
PYTHON_SCI ?= conda run -n duckfig python3
check-lr-binom:
	$(PYTHON_SCI) pipeline/_lr_binom.py

# ── Module-fingerprint cross-language check ────────────────────────────
# figures/figS03/_module_inputs.py's digest is recomputed independently in R by
# hh21_hh1716_nb_glmm.R; if the two ever disagree every figS03 consumer locks out with a
# false positive. Nothing else exercised that agreement. Needs Rscript + the digest package.
check-module-inputs:
	$(PYTHON_SCI) figures/figS03/_module_inputs.py

# ── Same-side duplicate reference vs e3 ────────────────────────────────
# d2l/d2n read references/same_side_duplicate_genotypes.tsv instead of e3; this
# recomputes the evidence from the committed e3 and raises on any disagreement.
# (The cluster runner FC_08656/run_regenerate_figure_data.sh, not distributed, runs it on the fresh e3.)
check-same-side-duplicates:
	$(PYTHON_SCI) pipeline/validate_same_side_duplicates.py

# ── f5 --from-cached-matrix obs loader ─────────────────────────────────
# A barcode repeated across libraries is legal; one repeated within a library raises.
check-locus-obs:
	$(PYTHON_SCI) scripts/test_locus_region_counts_obs.py

# Alias: `make check` → everything that gates a merge.
check: check-derived-data check-portable-paths check-cluster-python check-lr-binom check-module-inputs check-same-side-duplicates check-locus-obs

# ── Install pre-commit hook ────────────────────────────────────────────
# Idempotent — safe to run multiple times.  Registers .githooks/ as the
# hooks path for this repo so pre-commit fires the derived-data check
# (and any future hooks we add).
install-hooks:
	git config core.hooksPath .githooks
	chmod +x .githooks/pre-commit
	@echo "Hooks path set to .githooks/"

# ── Shell-script lint ─────────────────────────────────────────────────
# Requires shellcheck (brew install shellcheck).
lint-scripts:
	shellcheck scripts/*.sh .githooks/pre-commit
