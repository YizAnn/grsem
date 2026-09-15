source(file.path(Sys.getenv("GRSEM_PROJECT_ROOT"), "tools", "bootstrap.R"))
## Summarize completed Phase 2 raw simulation results without rerunning fits.
out <- output_dir("simulations/results_phase2")
recovery <- read.csv(file.path(out, "phase2_recovery_raw.csv"),
                     stringsAsFactors = FALSE)
roles <- read.csv(file.path(out, "phase2_role_raw.csv"),
                  stringsAsFactors = FALSE)

recovery_summary <- do.call(rbind, lapply(split(recovery,
  interaction(recovery$scenario, recovery$method, drop = TRUE)), function(z) {
    ok <- z$status != "failed"
    data.frame(
      scenario = z$scenario[1], method = z$method[1], replications = nrow(z),
      failure_rate = mean(!ok), group_tpr = mean(z$group_tpr[ok], na.rm = TRUE),
      group_fpr = mean(z$group_fpr[ok], na.rm = TRUE),
      parameter_tpr = mean(z$parameter_tpr[ok], na.rm = TRUE),
      parameter_fpr = mean(z$parameter_fpr[ok], na.rm = TRUE),
      loading_rmse = mean(z$loading_rmse[ok], na.rm = TRUE),
      runtime_seconds = mean(z$runtime_seconds[ok], na.rm = TRUE)
    )
  }))

role_summary <- do.call(rbind, lapply(split(roles, roles$scenario), function(z) {
  true_status <- if (identical(z$scenario[1], "role_mimic"))
    z$mimic_status else z$problem_factor_status
  all_failed <- z$reflective_status == "failed" &
    z$mimic_status == "failed" & z$problem_factor_status == "failed"
  any_failed <- z$reflective_status == "failed" |
    z$mimic_status == "failed" | z$problem_factor_status == "failed"
  data.frame(
    scenario = z$scenario[1], replications = nrow(z),
    true_model_failure_rate = mean(true_status == "failed"),
    any_model_failure_rate = mean(any_failed),
    all_models_failure_rate = mean(all_failed),
    role_recovery_unconditional = mean(z$role_recovered %in% TRUE),
    role_recovery_given_comparison = mean(z$role_recovered, na.rm = TRUE),
    comparisons_available = sum(!is.na(z$role_recovered))
  )
}))

write.csv(recovery_summary,
          file.path(out, "phase2_recovery_summary.csv"), row.names = FALSE)
write.csv(role_summary, file.path(out, "phase2_role_summary.csv"),
          row.names = FALSE)
print(recovery_summary)
print(role_summary)
