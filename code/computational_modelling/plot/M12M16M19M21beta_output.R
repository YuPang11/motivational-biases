

# -------- 1) Directory of the original Stan results (unchanged) --------
dir_hc1 <- "F:/PIT/model/EEGmodelHC/HC/results/stan"
file_fit_rda <- file.path(dir_hc1, "M12fits_obsData.Rdata")

# -------- 2) Directory where you want to save the exported data --------
out_dir <- "H:/eegplot/M4b"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------- 3) Load f2 and export to CSV --------
e <- new.env()
load(file_fit_rda, envir = e)   # f2 is required (dList is optional)
stopifnot(exists("f2", envir = e))

f2mat <- as.matrix(e$f2)

out_csv <- file.path(out_dir, "M12_f2_draws.csv")
write.csv(f2mat, out_csv, row.names = FALSE)

message("Saved to: ", out_csv)




# -------- 1) Directory of the original Stan results (unchanged) --------
dir_hc1 <- "F:/PIT/model/EEGmodelHC/HC/results/stan"
file_fit_rda <- file.path(dir_hc1, "M16fits_obsData.Rdata")

# -------- 2) Directory where you want to save the exported data --------
out_dir <- "F:/PIT/ISPSplot"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------- 3) Load f2 and export to CSV --------
e <- new.env()
load(file_fit_rda, envir = e)   # f2 is required (dList is optional)
stopifnot(exists("f2", envir = e))

f2mat <- as.matrix(e$f2)

out_csv <- file.path(out_dir, "M16_f2_drawshc.csv")
write.csv(f2mat, out_csv, row.names = FALSE)

message("Saved to: ", out_csv)



# -------- 1) Directory of the original Stan results (unchanged) --------
dir_hc1 <- "F:/PIT/model/EEGmodelHC/HC/results/stan"
file_fit_rda <- file.path(dir_hc1, "M19fits_obsData.Rdata")

# -------- 2) Directory where you want to save the exported data --------
out_dir <- "F:/PIT/ISPSplot"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------- 3) Load f2 and export to CSV --------
e <- new.env()
load(file_fit_rda, envir = e)   # f2 is required (dList is optional)
stopifnot(exists("f2", envir = e))

f2mat <- as.matrix(e$f2)

out_csv <- file.path(out_dir, "M19_f2_drawshc.csv")
write.csv(f2mat, out_csv, row.names = FALSE)

message("Saved to: ", out_csv)



# -------- 1) Directory of the original Stan results (unchanged) --------
dir_hc1 <- "F:/PIT/model/EEGmodelHC/OCD/results/stan"
file_fit_rda <- file.path(dir_hc1, "M12fits_obsData.Rdata")

# -------- 2) Directory where you want to save the exported data --------
out_dir <- "H:/eegplot/M4b"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------- 3) Load f2 and export to CSV --------
e <- new.env()
load(file_fit_rda, envir = e)   # f2 is required (dList is optional)
stopifnot(exists("f2", envir = e))

f2mat <- as.matrix(e$f2)

out_csv <- file.path(out_dir, "M12_f2_drawsocd.csv")
write.csv(f2mat, out_csv, row.names = FALSE)

message("Saved to: ", out_csv)





# -------- 1) Directory of the original Stan results (unchanged) --------
dir_hc1 <- "F:/PIT/model/EEGmodelHC/OCD/results/stan"
file_fit_rda <- file.path(dir_hc1, "M16fits_obsData.Rdata")

# -------- 2) Directory where you want to save the exported data --------
out_dir <- "F:/PIT/ISPSplot"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------- 3) Load f2 and export to CSV --------
e <- new.env()
load(file_fit_rda, envir = e)   # f2 is required (dList is optional)
stopifnot(exists("f2", envir = e))

f2mat <- as.matrix(e$f2)

out_csv <- file.path(out_dir, "M16_f2_drawsocd.csv")
write.csv(f2mat, out_csv, row.names = FALSE)

message("Saved to: ", out_csv)



# -------- 1) Directory of the original Stan results (unchanged) --------
dir_hc1 <- "F:/PIT/model/EEGmodelHC/OCD/results/stan"
file_fit_rda <- file.path(dir_hc1, "M21fits_obsData.Rdata")

# -------- 2) Directory where you want to save the exported data --------
out_dir <- "F:/PIT/ISPSplot"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------- 3) Load f2 and export to CSV --------
e <- new.env()
load(file_fit_rda, envir = e)   # f2 is required (dList is optional)
stopifnot(exists("f2", envir = e))

f2mat <- as.matrix(e$f2)

out_csv <- file.path(out_dir, "M21_f2_drawsocd.csv")
write.csv(f2mat, out_csv, row.names = FALSE)

message("Saved to: ", out_csv)
