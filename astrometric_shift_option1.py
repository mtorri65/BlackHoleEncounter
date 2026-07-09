# --- Build a per-visit dataframe for binning ---
dfv = pd.DataFrame({
    "timestamp_utc": sched["timestamp_utc"],
    "night_date_utc": sched["timestamp_utc"].dt.floor("D"),
    "dRA_mas": dRA_mas,
    "dDec_mas": dDec_mas,
    "shift_mas": shift_mas,
})

# Count visits per night (should be 2 for your paired schedule, but robust)
nightly = (
    dfv.groupby("night_date_utc", as_index=False)
       .agg(
           n_visits=("shift_mas", "size"),
           mean_shift_mas=("shift_mas", "mean"),
           mean_dRA_mas=("dRA_mas", "mean"),
           mean_dDec_mas=("dDec_mas", "mean"),
           max_shift_mas=("shift_mas", "max"),
       )
       .sort_values("night_date_utc")
       .reset_index(drop=True)
)

# Optional: effective nightly 3σ threshold if you average visits that night
# (sigma reduces ~ 1/sqrt(N_visits))
if np.isfinite(sigma_ast):
    nightly["sigma_eff_mas"] = sigma_ast / np.sqrt(nightly["n_visits"])
    nightly["det_line_3sig_eff_mas"] = K_SIGMA * nightly["sigma_eff_mas"]
