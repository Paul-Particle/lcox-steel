"""Physical constants shared across the pipeline.

These are fixed physical values — not assumptions or config knobs. Add new
constants here rather than inlining magic numbers in scripts.
"""

# Lower heating value of hydrogen (kWh per kg, LHV)
H2_LHV_KWH_PER_KG: float = 33.33

# Hours in a nameplate year: turns a per-year capacity or quote into a per-hour
# rate. The report levelises over the run's own calendar years instead.
HOURS_PER_YEAR: float = 8760.0
