suppressMessages({
  library(readxl)
  library(CHNOSZ)
})

#function to calculate ZC from chemical composition
lipid_ZC <- function(C_cc, H_cc, N_cc, O_cc, P_cc, Zplus_cc, Zminus_cc, S_cc, Z=NA){
  this_charge <- ifelse(is.na(Z), Zplus_cc - Zminus_cc, Z) # assumes Zminus_cc is given as a positive number
  this_ZC <- (-(H_cc) + (2 * O_cc) + (3 * N_cc) -
      (5 * P_cc) - (4 * S_cc) + this_charge) / (C_cc)
  return(this_ZC)
}

# function to output a chemical formula from elemental & charge composition
lipid_formula <- function(C_cc, H_cc, N_cc, O_cc, P_cc, Zplus_cc, Zminus_cc, S_cc){
  default_scipen_setting <- getOption("scipen")
  options(scipen = 999)
  formula <- paste(
    ifelse(C_cc > 0, paste("C", ifelse(C_cc != 1, C_cc, ""), sep=""), ""),
    ifelse(H_cc > 0, paste("H", ifelse(H_cc != 1, H_cc, ""), sep=""), ""),
    ifelse(N_cc > 0, paste("N", ifelse(N_cc != 1, N_cc, ""), sep=""), ""),
    ifelse(O_cc > 0, paste("O", ifelse(O_cc != 1, O_cc, ""), sep=""), ""),
    ifelse(P_cc > 0, paste("P", ifelse(P_cc != 1, P_cc, ""), sep=""), ""),
    ifelse(S_cc > 0, paste("S", ifelse(S_cc != 1, S_cc, ""), sep=""), ""),
    ifelse(Zplus_cc > 0, paste("+", ifelse(Zplus_cc != 1, Zplus_cc, ""), sep=""), ""),
    ifelse(Zminus_cc > 0, paste("-", ifelse(Zminus_cc != 1, Zminus_cc, ""), sep=""), ""),
    sep="")
  options(scipen = default_scipen_setting)
  return(formula)
}

# function to retrieve element abundance from a formula
# (e.g. number of carbons in "C6H12")
# input is string (e.g. "C6H12"), string (e.g. "C")
get_element_abund <- function(formula, element){
  if(element %in% names(makeup(formula))){
    makeup(formula)[element]
  } else {
    0
  }
}

# function to increase the capacity of is.nan() to work with data frames
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}


# function to report weighted formulas
ave_formula <- function(C_ave, H_ave, O_ave, N_ave, P_ave, S_ave, Z_ave){

  ave_formula <- paste0("C", signif(C_ave, 3),
                        "H", signif(H_ave, 3),
                        "O", signif(O_ave, 3),
                        "N", signif(N_ave, 3),
                        "P", signif(P_ave, 3),
                        "S", signif(S_ave, 3),
                        "Z", signif(Z_ave, 3))

  return(ave_formula)
}

# function to remove samples that do not have peak areas for a specified IPL.
# Used only during creation of "mergedIPL"
rmv_smpl <- function(x){
  if(exists("IPLareasRFmol")){
    abund <- IPLareasRFmol
  } else {
    abund <- mipl[["IPLareasRFmol"]]
  }
  if(is.matrix(x) | is.data.frame(x)){
    as.matrix(x[, colSums(abund) != 0, drop = FALSE])
  } else if (is.vector(x)){
    as.matrix(x[colSums(abund) != 0, drop = FALSE])
  } else {
    print(paste("Error removing samples without peak areas for IPL",
      this_IPL, subcategory, subgroup, "during creation of mergedIPL"))
      # note: variable this_IPL is defined inside of IPL_process()
  }
}

# Function to extract IPLs from Agilent MassHunter Quant ouput and process peak areas.
# Input data are in data/IPL_data.xlsx
# Output is an R list object called IPL_master
IPL_process <- function(IPL_grouping, workbook_name, quant_sample_order,
            IPLworkbook_directory, subgroup_colname, most_abund_IPL_only=FALSE) {
  time <- proc.time()
  IPLs_to_display <- unlist(IPL_grouping)
  IPL_master <- list()
  skipped_IPL <- c()

  # load response factor (RF) table from Excel workbook
  suppressMessages({
    RF_table <- as.data.frame(read_excel(paste0(IPLworkbook_directory, workbook_name), sheet = "RF"))
  })

  # load IPL headgroup information from Excel workbook
  suppressMessages({
    wkst_head <- as.data.frame(read_excel(paste0(IPLworkbook_directory, workbook_name), sheet = "IPL properties"))
  })

  for (n in 1:length(IPLs_to_display)) {

    this_IPL <- IPLs_to_display[n]

    subgroup <- ""
    if (grepl("sub:", this_IPL)){
      subgroup <- gsub(".+ sub: ", "", this_IPL)
      this_IPL <- gsub(" sub: .+", "", this_IPL)
    }

    subcategory <- ""
    if (grepl("subcat:", this_IPL)){
      subcategory <- gsub(".+ subcat: ", "", this_IPL)
      this_IPL <- gsub(" subcat: .+", "", this_IPL)
    }


    ############################ IMPORTING #################################
    # Import IPL sheet from Excel workbook
    suppressMessages({
      wkst <- as.data.frame(read_excel(paste0(IPLworkbook_directory, workbook_name), sheet = this_IPL))
    })

    if(is.null(wkst[1, subgroup_colname])){
      stop(paste("Error in subgroup: no column called", subgroup_colname, "in worksheet", this_IPL))
    }

    if(is.null(wkst[1, "Subcategory"]) & nchar(subcategory) > 0){
      stop(paste("Error in subcategory: no 'Subcategory' column in worksheet", this_IPL))
    }

    if (this_IPL != "AR-IPLs") {
      # Assign response factors for IPL batches (excluding AR-IPLs, which are handled later)
      RF_batch1 <- 1/as.numeric(RF_table[RF_table$IPL==this_IPL, "batch1"])
      #RF_batch2 <- 1/as.numeric(RF_table[RF_table$IPL==this_IPL, "batch2"])
    }

    ############################ GROUPING ##################################
    ### UPDATE ONLY IF SAMPLES ARE ADDED

    # Column order
    cols_batch1 <- 1:6 # Cols 1-14 are batch1
    #cols_batch2 <- 15:18 # Cols 15-18 are batch2

    ### Get only IPL names, m/z, RTs (retention times), peak areas, and MI (manual integration flag)
    IPL <- wkst[1:nrow(wkst), 1:20]

    ### Ignore ISTD in PC-DAG
    #if (this_IPL == "PC-DAG") {
      #IPL <- IPL[1:nrow(IPL)-1, ] # Remove the last row in PC-DAG, which is an ISTD
    #}

    ############################ PROCESSING ################################
    ### Assign name to MZ column
    colnames(IPL)[2] <- "MZ"

    ### Rename '...1' columns introduced by readxl v1.3.1 using the column name to the left of each
    for (j in 1:ncol(IPL)) {
      if (grepl("\\.\\.\\.[[:digit:]]+", colnames(IPL)[j])) {
        colnames(IPL)[j] <- colnames(IPL)[j - 1]
      }
    }

    ### Create a matrix 'IPLareas' that will store all peak areas
    # nrow is -1 due to the two lines of headers in the csv
    IPLareas <- matrix(nrow = nrow(IPL)-1, ncol=0)
    rownames(IPLareas) <- IPL[seq(2, length(IPL[, 1])), 1] # names the rows
    class(IPLareas) <- "numeric" # convert from string to numeric


    # Fill 'IPLareas' with peak area columns
    for (j in 1:ncol(IPL)) {
      if (grepl("Area",IPL[1, j])) { # regex "Area" in row1 header
        bindme <- IPL[seq(2, length(IPL[, j])), j]
        bindme <- as.numeric(bindme)
        IPLareas <- cbind(IPLareas, bindme)
        colnames(IPLareas)[ncol(IPLareas)] <- colnames(IPL)[j]
      }
    }

    # Convert all NAs to zeroes
    IPLareas[is.na(IPLareas)] <- 0

    # Rename IPLareas colnames
    colnames(IPLareas) <- quant_sample_order

    ### Create a matrix 'IPLmi' that stores manual-integration(MI) T/F flags
    IPLmi <- matrix(nrow = nrow(IPL) - 1, ncol = 0) # nrow(IPL)-1
    rownames(IPLmi) <- IPL[seq(2, length(IPL[, 1])), 1] # names the rows

    # Fill 'IPLmi' with MI flag columns
    for (j in 1:ncol(IPL)) {
      if (grepl("MI", IPL[1, j])) { # regex "MI" in row1 header
        bindme <- IPL[seq(2, length(IPL[, j])), j]
        IPLmi <- cbind(IPLmi, bindme)
        colnames(IPLmi)[ncol(IPLmi)] <- colnames(IPL)[j]
      }
    }

    # imports MI columns as "False" (e.g in 1G-DAG) or "FALSE"
    # (e.g. in AR-IPLs) for some reason. The following line sets either to 0:
    IPLareas[IPLmi == "False" | IPLmi == "FALSE"] <- 0

    ##### multiply peak areas by RF
    ### Create a new matrix 'IPLareasRF' to contain RF-modified peak areas
    IPLareasRF <- IPLareas
    IPL_RFith <- IPLareas/IPLareas
    IPL_RFith[is.nan(IPL_RFith)] <- 1

    if (this_IPL != "AR-IPLs") {
      # Cols 1-14 are batch1 and Cols 15-18 are from batch2
      IPLareasRF[, cols_batch1] <- IPLareasRF[, cols_batch1] * RF_batch1
      #IPLareasRF[, cols_batch2] <- IPLareasRF[, cols_batch2] * RF_batch2
      IPL_RFith[, cols_batch1] <- IPL_RFith[, cols_batch1] * RF_batch1
      #IPL_RFith[, cols_batch2] <- IPL_RFith[, cols_batch2] * RF_batch2
    } else {
      # Assign response factors to AR-IPLs
      RF_batch1 <- c()
      #RF_batch2 <- c()
      for (j in 1:nrow(IPLareasRF)) {
        IPLareasRF[j, cols_batch1] <- IPLareasRF[j, cols_batch1] /
          as.numeric(wkst[j + 1, "batch1"])
        #IPLareasRF[j, cols_batch2] <- IPLareasRF[j, cols_batch2] /
          #as.numeric(wkst[j + 1, "batch2"])
        IPL_RFith[j, cols_batch1] <- IPL_RFith[j, cols_batch1] /
          as.numeric(wkst[j + 1, "batch1"])
       # IPL_RFith[j, cols_batch2] <- IPL_RFith[j, cols_batch2] /
          #as.numeric(wkst[j + 1, "batch2"])
        RF_batch1 <- c(RF_batch1, 1 / as.numeric(wkst[j + 1, "batch1"]))
        #RF_batch2 <- c(RF_batch2, 1 / as.numeric(wkst[j + 1, "batch2"]))
      }
    }

    ### divide peak areas by molar mass to get from wt* to mole*

    # create numeric vector 'masses' from m/z column of the imported IPL csv
    masses <- as.numeric(IPL[seq(2, nrow(IPL)), 2])

    # divide peak wt* by IPL m/z to get mole* IPL
    IPLareasRFmol <- sweep(IPLareasRF, MARGIN = 1, masses, '/')

    # GB Oct 2, 2020: this provides a method to check only the most abundant lipids at each site.
    # Not currently used, but good for additional insight.
    if(most_abund_IPL_only){
      for(col in colnames(IPLareasRFmol)){
        IPLareasRFmol[, col] <- floor(IPLareasRFmol[, col] / max(IPLareasRFmol[, col]))
        IPLareasRFmol[is.nan(IPLareasRFmol)] <- 0
      }
    }

    ############################ ZC, nC, nUnsat ############################

    ### Create a vector of chain lengths
    nC <- c()

    if (this_IPL != "GDGTs" && this_IPL != "AR-IPLs") {
      for (j in 1:nrow(IPLareasRFmol)) {
        # get nC from row name
        nC_raw <- as.numeric(substr(rownames(IPLareasRFmol)[j], 2, 3))
        # add adjusted nC to vector
        nC <- c(nC, nC_raw)
      }
    } else if (this_IPL == "GDGTs") {
      for (j in 1:nrow(IPLareasRFmol)) {
        nC <- c(nC, 80) # 80 carbons within GDGT chains
      }
    } else if (this_IPL == "AR-IPLs") {
      for (j in 1:nrow(IPLareasRFmol)) {
        nC <- c(nC, wkst[j + 1, "nC"])
      }
    }

    # create an nC_raw vector used later for calculation of ave. chain nC
    nC_raw <- nC

    ### Create a matrix containing nC of the ith IPL
    IPL_nCith <- IPLareasRFmol
    IPL_nCith[, seq(1, ncol(IPL_nCith))] <- nC

    ### Initialize vectors of chain unsaturations and rings
    nUnsat <- c()
    nPentRing <- c()
    nHexRing <- c()

    if (this_IPL != "GDGTs" && this_IPL != "AR-IPLs") {
      for (j in 1:nrow(IPLareasRFmol)) {
        # get nUnsat from row name and add to vector
        nUnsat <- c(nUnsat, as.numeric(substr(rownames(IPLareasRFmol)[j], 5, 6)))
        nPentRing <- c(nPentRing, 0)
        nHexRing <- c(nHexRing, 0)
      }
    } else if (this_IPL == "GDGTs") {
      for (j in 1:nrow(IPLareasRFmol)) {
        nUnsat <- c(nUnsat, 0) # Assuming unsaturations != rings
        nPentRing <- c(nPentRing, wkst[j + 1, "internal pentacyclic rings"])
        nHexRing <- c(nHexRing, wkst[j + 1, "internal hexacyclic rings"])
      }
    } else if (this_IPL == "AR-IPLs") {
      for (j in 1:nrow(IPLareasRFmol)) {
        nUnsat <- c(nUnsat, wkst[j + 1, "unsat"])
        nPentRing <- c(nPentRing, 0)
        nHexRing <- c(nHexRing, 0)
      }
    }

    nUnsat_raw <- nUnsat

    ### Create a matrix containing nUnsat of the ith IPL
    IPL_nUnsatith <- IPLareasRFmol
    IPL_nUnsatith[, seq(1, ncol(IPL_nUnsatith))] <- nUnsat

    ### Create a matrix containing nPentRing of the ith IPL
    IPL_nPentRingith <- IPLareasRFmol
    IPL_nPentRingith[, seq(1, ncol(IPL_nPentRingith))] <- nPentRing

    ### Create a matrix containing nHexRing of the ith IPL
    IPL_nHexRingith <- IPLareasRFmol
    IPL_nHexRingith[, seq(1, ncol(IPL_nHexRingith))] <- nHexRing

    ### Calculate ZCs of structures based on elemental abundances

    # If an IPL has a user-specified subgroup, check for
    # currently only works for AR-IPLs and GDGTs using 'subgroup' column in excel workbook

    if (length(intersect(subcategory, wkst_head[, 1])) == 0 & nchar(subcategory) > 0) {
      stop(paste("Error in subcategory: chosen subcategory '", subcategory, "' does not appear in IPL properties worksheet", sep=""))
    }

    if (subcategory == "") {
      head_lookup_IPL <- paste("^", this_IPL, "$", sep="")
    } else if (nchar(subcategory) > 0) {
      head_lookup_IPL <- paste("^", subcategory, "$", sep="")
    } else {
      print("Error in subcategory")
    }

    head_row <- which(grepl(head_lookup_IPL, wkst_head[, 1]))

    # create a vector of n_chains and n_chain_types (ester, ether, amide, or nonlinkage)
    n_chains <- as.numeric(wkst_head[head_row, "n_chains"])
    n_chains_ether <- as.numeric(wkst_head[head_row, "n_ether_linkages"])
    n_chains_amide <- as.numeric(wkst_head[head_row, "n_amide_linkages"])
    n_chains_ester <- as.numeric(wkst_head[head_row, "n_ester_linkages"])
    n_chains_nonlinkage <- as.numeric(wkst_head[head_row, "n_nonlinkages"])
    n_chains_hydroxylated <- as.numeric(wkst_head[head_row, "n_chain_hydroxylations"])
    n_chains_GDGT <- as.numeric(wkst_head[head_row, "n_chain_GDGT"])
    n_chains_AR <- as.numeric(wkst_head[head_row, "n_chain_AR"])

    # create vectors of number of heads and backbones
    n_heads <- as.numeric(wkst_head[head_row, "n_heads"])
    n_backbones <- as.numeric(wkst_head[head_row, "n_backbones"])

    # for each row (IPL) in an IPL worksheet, calc. chemical composition and ZC
    cc <- "Chemical composition:"
    row_cc <- seq(1, length(nC)) + 1
    C_cc      <- as.numeric(wkst[row_cc, paste(cc, "C")])
    H_cc      <- as.numeric(wkst[row_cc, paste(cc, "H")])
    N_cc      <- as.numeric(wkst[row_cc, paste(cc, "N")])
    O_cc      <- as.numeric(wkst[row_cc, paste(cc, "O")])
    P_cc      <- as.numeric(wkst[row_cc, paste(cc, "P")])
    Zplus_cc  <- as.numeric(wkst[row_cc, paste(cc, "+")])
    Zminus_cc <- as.numeric(wkst[row_cc, paste(cc, "-")])
    S_cc      <- as.numeric(wkst[row_cc, paste(cc, "S")])

    ZCith <- lipid_ZC(C_cc, H_cc, N_cc, O_cc, P_cc, Zplus_cc, Zminus_cc, S_cc)
    ZCith_formula <- lipid_formula(C_cc, H_cc, N_cc, O_cc, P_cc, Zplus_cc, Zminus_cc, S_cc)

    #subtract out headgroups then calculate lipid ZC
    cc <- "Headgroup"
    C_head      <- as.numeric(wkst_head[head_row, paste(cc, "C")])
    H_head      <- as.numeric(wkst_head[head_row, paste(cc, "H")])
    N_head      <- as.numeric(wkst_head[head_row, paste(cc, "N")])
    O_head      <- as.numeric(wkst_head[head_row, paste(cc, "O")])
    P_head      <- as.numeric(wkst_head[head_row, paste(cc, "P")])
    Zplus_head  <- as.numeric(wkst_head[head_row, paste(cc, "+")])
    Zminus_head <- as.numeric(wkst_head[head_row, paste(cc, "-")])
    S_head      <- as.numeric(wkst_head[head_row, paste(cc, "S")])

    cc <- "Backbone"
    C_bb <- as.numeric(wkst_head[head_row, paste(cc, "C")])
    H_bb <- as.numeric(wkst_head[head_row, paste(cc, "H")])
    N_bb <- as.numeric(wkst_head[head_row, paste(cc, "N")])
    O_bb <- as.numeric(wkst_head[head_row, paste(cc, "O")])

    C_all_chains    <- C_cc - C_head - C_bb
    H_all_chains    <- H_cc - H_head - H_bb
    N_all_chains    <- N_cc - N_head - N_bb
    O_all_chains    <- O_cc - O_head - O_bb

    IPL_all_head_formula <- lipid_formula(C_head, H_head, N_head, O_head, P_head, Zplus_head, Zminus_head, S_head)
    IPL_all_backbone_formula <- lipid_formula(C_bb, H_bb, N_bb, O_bb, 0, 0, 0, 0)
    IPL_all_chain_formula <- lipid_formula(C_all_chains, H_all_chains, N_all_chains, O_all_chains, 0, 0, 0, 0)

    head_referenced <- wkst_head[head_row, 1]

    ### Create a matrix containing ZCs of the ith IPL
    IPL_ZCith <- IPLareasRFmol
    IPL_ZCith_formula <- IPLareasRFmol
    IPL_all_head_formulaith <- IPLareasRFmol
    IPL_all_backbone_formulaith <- IPLareasRFmol
    IPL_all_chain_formulaith <- IPLareasRFmol

    rows <- seq(1, nrow(IPL_ZCith))
    IPL_ZCith[rows, ] <- ZCith
    IPL_ZCith_formula[rows, ] <- ZCith_formula
    IPL_all_head_formulaith[rows, ] <- IPL_all_head_formula
    IPL_all_backbone_formulaith[rows, ] <- IPL_all_backbone_formula
    IPL_all_chain_formulaith[rows, ] <- IPL_all_chain_formula

    ##################### Subgroup and Subcategory-related processing ##########################
    # if a subgroup and/or subcategory is specified...
    if (nchar(subgroup) > 0 | nchar(subcategory) > 0) {
      # delete all rows from current IPL matrices that do not match subgroup or subcategory
      subgroup_rows <- c()
      subcategory_rows <- c()
      sub_subcat_rows <- c()

      for (nn in 1:nrow(IPLareasRFmol) + 1) {
        if (nchar(subgroup) == 0) {
          subgroup_rows <- c(subgroup_rows, nn)
        } else {
          if(wkst[nn, subgroup_colname] == subgroup) {
            subgroup_rows <- c(subgroup_rows, nn)
          }
        }

        if (nchar(subcategory) == 0) {
          subcategory_rows <- c(subcategory_rows, nn)
        } else {
          if (wkst[nn, "Subcategory"] == subcategory) {
            subcategory_rows <- c(subcategory_rows, nn)
          }
        }

      }

      sub_subcat_rows <- Reduce(intersect, list(subgroup_rows, subcategory_rows))

      sub_subcat_rows <- sub_subcat_rows - 1

      IPLmi <- IPLmi[sub_subcat_rows, , drop = FALSE]
      IPLareas <- IPLareas[sub_subcat_rows, , drop = FALSE]
      masses <- masses[sub_subcat_rows]
      IPLareasRF <- IPLareasRF[sub_subcat_rows, , drop = FALSE]
      IPLareasRFmol <- IPLareasRFmol[sub_subcat_rows, , drop = FALSE]
      nC <- nC[sub_subcat_rows]
      nC_raw <- nC_raw[sub_subcat_rows]
      IPL_nCith <- IPL_nCith[sub_subcat_rows, , drop = FALSE]
      nUnsat <- nUnsat[sub_subcat_rows]
      nUnsat_raw <- nUnsat_raw[sub_subcat_rows]
      nPentRing <- nPentRing[sub_subcat_rows]
      nHexRing <- nHexRing[sub_subcat_rows]
      IPL_nUnsatith <- IPL_nUnsatith[sub_subcat_rows, , drop = FALSE]
      IPL_nPentRingith <- IPL_nPentRingith[sub_subcat_rows, , drop = FALSE]
      IPL_nHexRingith <- IPL_nHexRingith[sub_subcat_rows, , drop = FALSE]
      ZCith <- ZCith[sub_subcat_rows]
      IPL_ZCith <- IPL_ZCith[sub_subcat_rows, , drop = FALSE]
      IPL_ZCith_formula <- IPL_ZCith_formula[sub_subcat_rows, , drop = FALSE]

      IPL_all_head_formulaith <- IPL_all_head_formulaith[sub_subcat_rows, , drop = FALSE]
      IPL_all_backbone_formulaith <- IPL_all_backbone_formulaith[sub_subcat_rows, , drop = FALSE]
      IPL_all_chain_formulaith <- IPL_all_chain_formulaith[sub_subcat_rows, , drop = FALSE]

      RF_batch1 <- RF_batch1[sub_subcat_rows]
      #RF_batch2 <- RF_batch2[sub_subcat_rows]
      IPL_RFith <- IPL_RFith[sub_subcat_rows, , drop = FALSE]

      if (nchar(subgroup) > 0){
        subgroup <- paste(" sub: ", subgroup, sep = "")
      }

      if (nchar(subcategory) > 0){
        subcategory <- paste(" subcat: ", subcategory, sep = "")
      }
    }


    # If an IPLs has had all of its rows deleted (no subcategory or subgroup match) then
    # skip the rest of the processing and move onto the next IPL
    if (length(IPLareasRFmol) == 0 & length(IPLs_to_display) > 1){
      skipped_IPL <- c(skipped_IPL, paste(this_IPL, subcategory, subgroup, sep=""))
      next
    } else if (length(IPLareasRFmol) == 0 & length(IPLs_to_display) == 1){
      stop("No IPLs match chosen subgroup")
    }

    ################## Write processed data to IPLmaster #######################
    IPL_name <- paste(this_IPL, subcategory, subgroup, sep = "")

    IPL_master[[IPL_name]][["IPLmi"]] <- IPLmi
    IPL_master[[IPL_name]][["IPLareas"]] <- IPLareas
    IPL_master[[IPL_name]][["masses"]] <- masses
    IPL_master[[IPL_name]][["IPLareasRF"]] <- IPLareasRF
    IPL_master[[IPL_name]][["IPLareasRFmol"]] <- IPLareasRFmol
    IPL_master[[IPL_name]][["nC"]] <- nC
    IPL_master[[IPL_name]][["nC_raw"]] <- nC_raw
    IPL_master[[IPL_name]][["IPL_nCith"]] <- IPL_nCith
    IPL_master[[IPL_name]][["nUnsat"]] <- nUnsat
    IPL_master[[IPL_name]][["nUnsat_raw"]] <- nUnsat_raw
    IPL_master[[IPL_name]][["nPentRing"]] <- nPentRing
    IPL_master[[IPL_name]][["nHexRing"]] <- nHexRing
    IPL_master[[IPL_name]][["IPL_nUnsatith"]] <- IPL_nUnsatith
    IPL_master[[IPL_name]][["IPL_nPentRingith"]] <- IPL_nPentRingith
    IPL_master[[IPL_name]][["IPL_nHexRingith"]] <- IPL_nHexRingith
    IPL_master[[IPL_name]][["ZCith"]] <- ZCith
    IPL_master[[IPL_name]][["IPL_ZCith"]] <- IPL_ZCith
    IPL_master[[IPL_name]][["IPL_ZCith_formula"]] <- IPL_ZCith_formula
    IPL_master[[IPL_name]][["IPL_all_head_formulaith"]] <- IPL_all_head_formulaith
    IPL_master[[IPL_name]][["IPL_all_backbone_formulaith"]] <- IPL_all_backbone_formulaith
    IPL_master[[IPL_name]][["IPL_all_chain_formulaith"]] <- IPL_all_chain_formulaith
    IPL_master[[IPL_name]][["n_chains"]] <- n_chains
    IPL_master[[IPL_name]][["n_chains_ether"]] <- n_chains_ether
    IPL_master[[IPL_name]][["n_chains_amide"]] <- n_chains_amide
    IPL_master[[IPL_name]][["n_chains_ester"]] <- n_chains_ester
    IPL_master[[IPL_name]][["n_chains_nonlinkage"]] <- n_chains_nonlinkage
    IPL_master[[IPL_name]][["n_chains_hydroxylated"]] <- n_chains_hydroxylated
    IPL_master[[IPL_name]][["n_chains_GDGT"]] <- n_chains_GDGT
    IPL_master[[IPL_name]][["n_chains_AR"]] <- n_chains_AR
    IPL_master[[IPL_name]][["IPL_RFith"]] <- IPL_RFith

    ################# Abundance-weighted ZC, nC, nUnsat ########################
    ### June 28th - re-evaluation of per-chain calculations
    # create a matrix of nchains (number of chains in the ith IPL)
    n_chains_matrix <- IPL_nCith
    n_chains_ether_matrix <- IPL_nCith
    n_chains_ester_matrix <- IPL_nCith
    n_chains_amide_matrix <- IPL_nCith
    n_chains_nonlinkage_matrix <- IPL_nCith
    n_chains_hydroxylated_matrix <- IPL_nCith
    n_chains_GDGT_matrix <- IPL_nCith
    n_chains_AR_matrix <- IPL_nCith
    n_heads_matrix <- IPL_nCith
    n_backbones_matrix <- IPL_nCith
    n_chains_matrix[seq(1,nrow(n_chains_matrix)),seq(1,ncol(n_chains_matrix))] <- n_chains
    n_chains_ether_matrix[seq(1,nrow(n_chains_matrix)),seq(1,ncol(n_chains_matrix))] <- n_chains_ether
    n_chains_ester_matrix[seq(1,nrow(n_chains_matrix)),seq(1,ncol(n_chains_matrix))] <- n_chains_ester
    n_chains_amide_matrix[seq(1,nrow(n_chains_matrix)),seq(1,ncol(n_chains_matrix))] <- n_chains_amide
    n_chains_nonlinkage_matrix[seq(1,nrow(n_chains_matrix)),seq(1,ncol(n_chains_matrix))] <- n_chains_nonlinkage
    n_chains_hydroxylated_matrix[seq(1,nrow(n_chains_matrix)),seq(1,ncol(n_chains_matrix))] <- n_chains_hydroxylated
    n_chains_GDGT_matrix[seq(1,nrow(n_chains_matrix)),seq(1,ncol(n_chains_matrix))] <- n_chains_GDGT
    n_chains_AR_matrix[seq(1,nrow(n_chains_matrix)),seq(1,ncol(n_chains_matrix))] <- n_chains_AR
    n_heads_matrix[seq(1,nrow(n_heads_matrix)),seq(1,ncol(n_heads_matrix))] <- n_heads
    n_backbones_matrix[seq(1,nrow(n_backbones_matrix)),seq(1,ncol(n_backbones_matrix))] <- n_backbones
    IPL_master[[IPL_name]][["n_chains_matrix"]] <- n_chains_matrix
    IPL_master[[IPL_name]][["n_chains_ether_matrix"]] <- n_chains_ether_matrix
    IPL_master[[IPL_name]][["n_chains_ester_matrix"]] <- n_chains_ester_matrix
    IPL_master[[IPL_name]][["n_chains_amide_matrix"]] <- n_chains_amide_matrix
    IPL_master[[IPL_name]][["n_chains_nonlinkage_matrix"]] <- n_chains_nonlinkage_matrix
    IPL_master[[IPL_name]][["n_chains_hydroxylated_matrix"]] <- n_chains_hydroxylated_matrix
    IPL_master[[IPL_name]][["n_chains_GDGT_matrix"]] <- n_chains_GDGT_matrix
    IPL_master[[IPL_name]][["n_chains_AR_matrix"]] <- n_chains_AR_matrix
    IPL_master[[IPL_name]][["n_heads_matrix"]] <- n_heads_matrix
    IPL_master[[IPL_name]][["n_backbones_matrix"]] <- n_backbones_matrix



    ### Get abundance-weighted average properties of all IPLs at each site
    # nC
    IPL_nCave <- colSums(sweep(IPLareasRFmol, MARGIN = 1, nC_raw, "*",
      check.margin = FALSE)) / colSums(IPLareasRFmol * n_chains_matrix)

    # nUnsat
    IPL_nUnsatave <- colSums(sweep(IPLareasRFmol, MARGIN = 1, nUnsat_raw, "*",
      check.margin = FALSE)) / colSums(IPLareasRFmol * n_chains_matrix)

    # nPentRing
    IPL_nPentRingave <- colSums(sweep(IPLareasRFmol, MARGIN = 1, nPentRing, "*",
      check.margin = FALSE)) / colSums(IPLareasRFmol * n_chains_matrix)

    # nHexRing
    IPL_nHexRingave <- colSums(sweep(IPLareasRFmol, MARGIN = 1, nHexRing, "*",
      check.margin = FALSE)) / colSums(IPLareasRFmol * n_chains_matrix)

    # ZC
    IPL_ZCave <- colSums(sweep(IPLareasRFmol, MARGIN = 1, ZCith, "*", check.margin = FALSE)) / colSums(IPLareasRFmol)

    # replace NaNs with 0
    IPL_nCave[is.nan(IPL_nCave)] <- 0
    IPL_nUnsatave[is.nan(IPL_nUnsatave)] <- 0
    IPL_nPentRingave[is.nan(IPL_nPentRingave)] <- 0
    IPL_nHexRingave[is.nan(IPL_nHexRingave)] <- 0
    IPL_ZCave[is.nan(IPL_ZCave)] <- 0

    # Report data
    IPL_master[[IPL_name]][["IPL_nCave"]] <- IPL_nCave
    IPL_master[[IPL_name]][["IPL_nUnsatave"]] <- IPL_nUnsatave
    IPL_master[[IPL_name]][["IPL_nPentRingave"]] <- IPL_nPentRingave
    IPL_master[[IPL_name]][["IPL_nHexRingave"]] <- IPL_nHexRingave
    IPL_master[[IPL_name]][["IPL_ZCave"]] <- IPL_ZCave


    ################# Mean ZC, nC, nUnsat for plotting ########################
    # Remove samples with no peak areas from both x, y, and mean-y values
    y_ZCith <- as.matrix(IPL_ZCith[, colSums(IPLareasRFmol) != 0, drop = FALSE])
    y_mean_ZCith <- as.matrix(IPL_ZCave[colSums(IPLareasRFmol) != 0, drop = FALSE])

    y_nCith <-  as.matrix(IPL_nCith[, colSums(IPLareasRFmol) != 0, drop = FALSE])
    y_mean_nCith <- as.matrix(IPL_nCave[colSums(IPLareasRFmol) != 0, drop = FALSE])

    y_nUnsatith <-  as.matrix(IPL_nUnsatith[, colSums(IPLareasRFmol) != 0, drop = FALSE])
    y_mean_nUnsatith <- as.matrix(IPL_nUnsatave[colSums(IPLareasRFmol) != 0, drop = FALSE])

    y_nPentRingith <-  as.matrix(IPL_nPentRingith[, colSums(IPLareasRFmol) != 0, drop = FALSE])
    y_mean_nPentRingith <- as.matrix(IPL_nPentRingave[colSums(IPLareasRFmol) != 0, drop = FALSE])

    y_nHexRingith <-  as.matrix(IPL_nHexRingith[, colSums(IPLareasRFmol) != 0, drop = FALSE])
    y_mean_nHexRingith <- as.matrix(IPL_nHexRingave[colSums(IPLareasRFmol) != 0, drop = FALSE])

    wt_matrix <- as.matrix(IPLareasRFmol[, colSums(IPLareasRFmol) != 0, drop = FALSE])
    wt_for_mean <- as.matrix(colSums(wt_matrix))

    IPL_master[[IPL_name]][["y_ZCith"]] <- y_ZCith

    IPL_master[[IPL_name]][["y_nCith"]] <- y_nCith
    IPL_master[[IPL_name]][["y_nUnsatith"]] <- y_nUnsatith
    IPL_master[[IPL_name]][["y_nPentRingith"]] <- y_nPentRingith
    IPL_master[[IPL_name]][["y_nHexRingith"]] <- y_nHexRingith

    IPL_master[[IPL_name]][["y_mean_ZCith"]] <- y_mean_ZCith

    IPL_master[[IPL_name]][["y_mean_nCith"]] <- y_mean_nCith
    IPL_master[[IPL_name]][["y_mean_nUnsatith"]] <- y_mean_nUnsatith
    IPL_master[[IPL_name]][["y_mean_nPentRingith"]] <- y_mean_nPentRingith
    IPL_master[[IPL_name]][["y_mean_nHexRingith"]] <- y_mean_nHexRingith

    IPL_master[[IPL_name]][["wt_matrix"]] <- wt_matrix
    IPL_master[[IPL_name]][["wt_for_mean"]] <- wt_for_mean

    # Store information about individual IPLs
    IPL_master[[IPL_name]][["info"]][["display"]] <- IPLs_to_display[n]
    IPL_master[[IPL_name]][["info"]][["RF"]] <- list(c(1/RF_batch1)) #,1/RF_batch2))
  }

  ################################################################
  ####################### CREATE mergedIPL #######################
  ################################################################

  # If there is grouping/merging, perform the merge:
  for (i in 1:length(IPL_grouping)) {

    mipl <- list()

    for (n in 1:length(IPL_grouping[[i]])) {
      this_IPL <- IPL_grouping[[i]][n]
      IPL_name <- IPL_grouping[[i]][n]

      mipl[["n_chains_matrix"]] <- rbind(mipl[["n_chains_matrix"]],
        IPL_master[[IPL_name]][["n_chains_matrix"]])
      mipl[["n_chains_ether_matrix"]] <- rbind(mipl[["n_chains_ether_matrix"]],
        IPL_master[[IPL_name]][["n_chains_ether_matrix"]])
      mipl[["n_chains_ester_matrix"]] <- rbind(mipl[["n_chains_ester_matrix"]],
        IPL_master[[IPL_name]][["n_chains_ester_matrix"]])
      mipl[["n_chains_amide_matrix"]] <- rbind(mipl[["n_chains_amide_matrix"]],
        IPL_master[[IPL_name]][["n_chains_amide_matrix"]])
      mipl[["n_chains_nonlinkage_matrix"]] <- rbind(mipl[["n_chains_nonlinkage_matrix"]],
        IPL_master[[IPL_name]][["n_chains_nonlinkage_matrix"]])
      mipl[["n_chains_hydroxylated_matrix"]] <- rbind(mipl[["n_chains_hydroxylated_matrix"]],
        IPL_master[[IPL_name]][["n_chains_hydroxylated_matrix"]])
      mipl[["n_chains_GDGT_matrix"]] <- rbind(mipl[["n_chains_GDGT_matrix"]],
        IPL_master[[IPL_name]][["n_chains_GDGT_matrix"]])
      mipl[["n_chains_AR_matrix"]] <- rbind(mipl[["n_chains_AR_matrix"]],
        IPL_master[[IPL_name]][["n_chains_AR_matrix"]])
      mipl[["n_heads_matrix"]] <- rbind(mipl[["n_heads_matrix"]],
        IPL_master[[IPL_name]][["n_heads_matrix"]])
      mipl[["n_backbones_matrix"]] <- rbind(mipl[["n_backbones_matrix"]],
        IPL_master[[IPL_name]][["n_backbones_matrix"]])

      mipl[["IPLareas"]] <- rbind(mipl[["IPLareas"]],
        IPL_master[[IPL_name]][["IPLareas"]])
      mipl[["IPLareasRFmol"]] <- rbind(mipl[["IPLareasRFmol"]],
        IPL_master[[IPL_name]][["IPLareasRFmol"]])
      mipl[["IPL_nCith"]] <- rbind(mipl[["IPL_nCith"]],
        IPL_master[[IPL_name]][["IPL_nCith"]])
      mipl[["IPL_nUnsatith"]] <- rbind(mipl[["IPL_nUnsatith"]],
        IPL_master[[IPL_name]][["IPL_nUnsatith"]])
      mipl[["IPL_nPentRingith"]] <- rbind(mipl[["IPL_nPentRingith"]],
        IPL_master[[IPL_name]][["IPL_nPentRingith"]])
      mipl[["IPL_nHexRingith"]] <- rbind(mipl[["IPL_nHexRingith"]],
        IPL_master[[IPL_name]][["IPL_nHexRingith"]])
      mipl[["IPL_ZCith"]] <- rbind(mipl[["IPL_ZCith"]],
        IPL_master[[IPL_name]][["IPL_ZCith"]])
      mipl[["ZCith"]] <- c(mipl[["ZCith"]],
        IPL_master[[IPL_name]][["ZCith"]])
      mipl[["nC"]] <- c(mipl[["nC"]],
        IPL_master[[IPL_name]][["nC"]])
      mipl[["nC_raw"]] <- c(mipl[["nC_raw"]],
        IPL_master[[IPL_name]][["nC_raw"]])
      mipl[["nUnsat"]] <- c(mipl[["nUnsat"]],
        IPL_master[[IPL_name]][["nUnsat"]])
        mipl[["nUnsat_raw"]] <- c(mipl[["nUnsat_raw"]],
          IPL_master[[IPL_name]][["nUnsat_raw"]])
      mipl[["nPentRing"]] <- c(mipl[["nPentRing"]],
        IPL_master[[IPL_name]][["nPentRing"]])
      mipl[["nHexRing"]] <- c(mipl[["nHexRing"]],
        IPL_master[[IPL_name]][["nHexRing"]])
      mipl[["masses"]] <- c(mipl[["masses"]],
        IPL_master[[IPL_name]][["masses"]])

      mipl[["n_chains"]] <- c(mipl[["n_chains"]],
        IPL_master[[IPL_name]][["n_chains"]])
      mipl[["n_chains_ether"]] <- c(mipl[["n_chains_ether"]],
        IPL_master[[IPL_name]][["n_chains_ether"]])
      mipl[["n_chains_ester"]] <- c(mipl[["n_chains_ester"]],
        IPL_master[[IPL_name]][["n_chains_ester"]])
      mipl[["n_chains_amide"]] <- c(mipl[["n_chains_amide"]],
        IPL_master[[IPL_name]][["n_chains_amide"]])
      mipl[["n_chains_nonlinkage"]] <- c(mipl[["n_chains_nonlinkage"]],
        IPL_master[[IPL_name]][["n_chains_nonlinkage"]])
      mipl[["n_chains_hydroxylated"]] <- c(mipl[["n_chains_hydroxylated"]],
        IPL_master[[IPL_name]][["n_chains_hydroxylated"]])
      mipl[["n_chains_GDGT"]] <- c(mipl[["n_chains_GDGT"]],
        IPL_master[[IPL_name]][["n_chains_GDGT"]])
      mipl[["n_chains_AR"]] <- c(mipl[["n_chains_AR"]],
        IPL_master[[IPL_name]][["n_chains_AR"]])

      mipl[["IPL_ZCith_formula"]] <- rbind(mipl[["IPL_ZCith_formula"]],
        IPL_master[[IPL_name]][["IPL_ZCith_formula"]])

      mipl[["IPL_all_head_formulaith"]] <- rbind(mipl[["IPL_all_head_formulaith"]],
        IPL_master[[IPL_name]][["IPL_all_head_formulaith"]])
      mipl[["IPL_all_backbone_formulaith"]] <- rbind(mipl[["IPL_all_backbone_formulaith"]],
        IPL_master[[IPL_name]][["IPL_all_backbone_formulaith"]])
      mipl[["IPL_all_chain_formulaith"]] <- rbind(mipl[["IPL_all_chain_formulaith"]],
        IPL_master[[IPL_name]][["IPL_all_chain_formulaith"]])

      mipl[["this IPL"]] <- c(mipl[["this IPL"]], IPL_name)
      mipl[["IPL_RFith nrow"]] <- c(mipl[["IPL_RFith nrow"]],
        nrow(IPL_master[[IPL_name]][["IPL_RFith"]]))
      mipl[["IPL_RFith"]] <- rbind(mipl[["IPL_RFith"]],
        IPL_master[[IPL_name]][["IPL_RFith"]])
    }

    ### Get abundance-weighted average nC of all IPLs at each site
    IPLareasRFmol <- mipl[["IPLareasRFmol"]]

    IPL_nCave_matrix <- IPLareasRFmol
    IPL_nUnsatave_matrix <- IPLareasRFmol
    IPL_nPentRingave_matrix <- IPLareasRFmol
    IPL_nHexRingave_matrix <- IPLareasRFmol
    IPL_ZCave_matrix <- IPLareasRFmol

    ZCith <- mipl[["ZCith"]]
    nC <- mipl[["nC"]]
    nC_raw <- mipl[["nC_raw"]]
    nUnsat <- mipl[["nUnsat"]]
    nUnsat_raw <- mipl[["nUnsat_raw"]]
    nPentRing <- mipl[["nPentRing"]]
    nHexRing <- mipl[["nHexRing"]]
    IPL_nCith <- mipl[["IPL_nCith"]]
    IPL_nUnsatith <- mipl[["IPL_nUnsatith"]]
    IPL_nPentRingith <- mipl[["IPL_nPentRingith"]]
    IPL_nHexRingith <- mipl[["IPL_nHexRingith"]]
    IPL_ZCith <- mipl[["IPL_ZCith"]]
    n_chains_matrix <- mipl[["n_chains_matrix"]]
    n_chains_ether_matrix <- mipl[["n_chains_ether_matrix"]]
    n_chains_ester_matrix <- mipl[["n_chains_ester_matrix"]]
    n_chains_amide_matrix <- mipl[["n_chains_amide_matrix"]]
    n_chains_hydroxylated_matrix <- mipl[["n_chains_hydroxylated_matrix"]]
    n_chains_nonlinkage_matrix <- mipl[["n_chains_nonlinkage_matrix"]]
    n_chains_GDGT_matrix <- mipl[["n_chains_GDGT_matrix"]]
    n_chains_AR_matrix <- mipl[["n_chains_AR_matrix"]]
    n_heads_matrix <- mipl[["n_heads_matrix"]]
    n_backbones_matrix <- mipl[["n_backbones_matrix"]]

    IPL_all_IPL_formulaith <- mipl[["IPL_ZCith_formula"]]
    IPL_all_head_formulaith <- mipl[["IPL_all_head_formulaith"]]
    IPL_all_backbone_formulaith <- mipl[["IPL_all_backbone_formulaith"]]
    IPL_all_chain_formulaith <- mipl[["IPL_all_chain_formulaith"]]

    # get formulas of components
    IPL_formulas <- IPL_all_IPL_formulaith[, 1]
    head_formulas <- IPL_all_head_formulaith[, 1]
    backbone_formulas <- IPL_all_backbone_formulaith[, 1]
    chain_formulas <- IPL_all_chain_formulaith[, 1]

    # # create matrices storing abundance of each element in each component of the ith IPL
    C_all_IPL <- sapply(IPL_formulas, FUN = get_element_abund, element = "C")
    H_all_IPL <- sapply(IPL_formulas, FUN = get_element_abund, element = "H")
    N_all_IPL <- sapply(IPL_formulas, FUN = get_element_abund, element = "N")
    O_all_IPL <- sapply(IPL_formulas, FUN = get_element_abund, element = "O")
    P_all_IPL <- sapply(IPL_formulas, FUN = get_element_abund, element = "P")
    S_all_IPL <- sapply(IPL_formulas, FUN = get_element_abund, element = "S")
    Z_all_IPL <- sapply(IPL_formulas, FUN = get_element_abund, element = "Z")
    C_all_head <- sapply(head_formulas, FUN = get_element_abund, element = "C")
    H_all_head <- sapply(head_formulas, FUN = get_element_abund, element = "H")
    N_all_head <- sapply(head_formulas, FUN = get_element_abund, element = "N")
    O_all_head <- sapply(head_formulas, FUN = get_element_abund, element = "O")
    P_all_head <- sapply(head_formulas, FUN = get_element_abund, element = "P")
    S_all_head <- sapply(head_formulas, FUN = get_element_abund, element = "S")
    Z_all_head <- sapply(head_formulas, FUN = get_element_abund, element = "Z")
    C_all_backbone <- sapply(backbone_formulas, FUN = get_element_abund, element = "C")
    H_all_backbone <- sapply(backbone_formulas, FUN = get_element_abund, element = "H")
    N_all_backbone <- sapply(backbone_formulas, FUN = get_element_abund, element = "N")
    O_all_backbone <- sapply(backbone_formulas, FUN = get_element_abund, element = "O")
    P_all_backbone <- sapply(backbone_formulas, FUN = get_element_abund, element = "P")
    S_all_backbone <- sapply(backbone_formulas, FUN = get_element_abund, element = "S")
    Z_all_backbone <- sapply(backbone_formulas, FUN = get_element_abund, element = "Z")
    C_all_chain <- sapply(chain_formulas, FUN = get_element_abund, element = "C")
    H_all_chain <- sapply(chain_formulas, FUN = get_element_abund, element = "H")
    N_all_chain <- sapply(chain_formulas, FUN = get_element_abund, element = "N")
    O_all_chain <- sapply(chain_formulas, FUN = get_element_abund, element = "O")
    P_all_chain <- sapply(chain_formulas, FUN = get_element_abund, element = "P")
    S_all_chain <- sapply(chain_formulas, FUN = get_element_abund, element = "S")
    Z_all_chain <- sapply(chain_formulas, FUN = get_element_abund, element = "Z")


    elements <- c("C", "H", "N", "O", "P", "S", "Z")
    for (element in elements){

      # create matrices from vectors
      eval(parse(text = paste0(element, "_all_IPL_matrix <- n_heads_matrix"))) # n_heads_matrix is okay here
      eval(parse(text = paste0(element, "_all_head_matrix <- n_heads_matrix")))
      eval(parse(text = paste0(element, "_all_backbone_matrix <- n_backbones_matrix")))
      eval(parse(text = paste0(element, "_all_chain_matrix <- n_chains_matrix")))

      eval(parse(text = paste0(element, "_all_IPL_matrix[seq(1, nrow(n_heads_matrix)), seq(1, ncol(n_heads_matrix))] <- ", element, "_all_IPL"))) # using n_heads matrix should be okay here
      eval(parse(text = paste0(element, "_all_head_matrix[seq(1, nrow(n_heads_matrix)), seq(1, ncol(n_heads_matrix))] <- ", element, "_all_head")))
      eval(parse(text = paste0(element, "_all_backbone_matrix[seq(1, nrow(n_backbones_matrix)), seq(1, ncol(n_backbones_matrix))] <- ", element, "_all_backbone")))
      eval(parse(text = paste0(element, "_all_chain_matrix[seq(1, nrow(n_chains_matrix)), seq(1, ncol(n_chains_matrix))] <- ", element, "_all_chain")))

      # Create a matrix containing per-component ith element for this IPL (not an average)
      eval(parse(text = paste0("IPL_", element, "_head <- ", element, "_all_head_matrix / n_heads_matrix")))
      eval(parse(text = paste0("IPL_", element, "_backbone <- ", element, "_all_backbone_matrix / n_backbones_matrix")))
      eval(parse(text = paste0("IPL_", element, "_chain <- ", element, "_all_chain_matrix / n_chains_matrix")))

      # replace nan with 0 in element_all_component_matrix
      eval(parse(text = paste0(element, "_all_IPL_matrix[is.nan(", element, "_all_IPL_matrix)] <- 0")))
      eval(parse(text = paste0(element, "_all_head_matrix[is.nan(", element, "_all_head_matrix)] <- 0")))
      eval(parse(text = paste0(element, "_all_backbone_matrix[is.nan(", element, "_all_backbone_matrix)] <- 0")))
      eval(parse(text = paste0(element, "_all_chain_matrix[is.nan(", element, "_all_chain_matrix)] <- 0")))

      # Calculate component abundance-weighted formula value of current element in loop.
      # Returns a vector labeled with sample names and value of current element
      eval(parse(text = paste0("IPL_", element, "_IPL_ave <- colSums(", element, "_all_IPL_matrix * IPLareasRFmol) / colSums(IPLareasRFmol)")))
      eval(parse(text = paste0("IPL_", element, "_head_ave <- colSums(", element, "_all_head_matrix * IPLareasRFmol) / colSums(n_heads_matrix * IPLareasRFmol)")))
      eval(parse(text = paste0("IPL_", element, "_backbone_ave <- colSums(", element, "_all_backbone_matrix * IPLareasRFmol) / colSums(n_backbones_matrix * IPLareasRFmol)")))
      eval(parse(text = paste0("IPL_", element, "_chain_ave <- colSums(", element, "_all_chain_matrix * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)")))

      # replace nan with 0 in IPL_element_component_ave matrix
      eval(parse(text = paste0("IPL_", element, "_IPL_ave[is.nan(IPL_", element, "_IPL_ave)] <- 0")))
      eval(parse(text = paste0("IPL_", element, "_head_ave[is.nan(IPL_", element, "_head_ave)] <- 0")))
      eval(parse(text = paste0("IPL_", element, "_backbone_ave[is.nan(IPL_", element, "_backbone_ave)] <- 0")))
      eval(parse(text = paste0("IPL_", element, "_chain_ave[is.nan(IPL_", element, "_chain_ave)] <- 0")))

    }




    # Create a matrix of IPL_<component>_formulaith, which differs from IPL_all_<component>_formulaith
    # in that IPL_<component>_formulaith is PER CHAIN and IPL_all_<component>_formulaith for ALL <component> OF IPL.
    # This matrix will be used for plotting.
    IPL_head_formulaith <- as.data.frame(IPL_C_head) # initialize matrix
    IPL_backbone_formulaith <- as.data.frame(IPL_C_backbone) # initialize matrix
    IPL_chain_formulaith <- as.data.frame(IPL_C_chain) # initialize matrix


    for (row in 1:nrow(IPL_C_chain)){
      IPL_head_formulaith[row, 1] <- paste0("C", IPL_C_head[row, 1], "H", IPL_H_head[row, 1], "O", IPL_O_head[row, 1], "N", IPL_N_head[row, 1], "P", IPL_P_head[row, 1], "S", IPL_S_head[row, 1], "Z", IPL_Z_head[row, 1])
      IPL_backbone_formulaith[row, 1] <- paste0("C", IPL_C_backbone[row, 1], "H", IPL_H_backbone[row, 1], "O", IPL_O_backbone[row, 1], "N", IPL_N_backbone[row, 1], "P", IPL_P_backbone[row, 1], "S", IPL_S_backbone[row, 1], "Z", IPL_Z_backbone[row, 1])
      IPL_chain_formulaith[row, 1] <- paste0("C", IPL_C_chain[row, 1], "H", IPL_H_chain[row, 1], "O", IPL_O_chain[row, 1], "N", IPL_N_chain[row, 1], "P", IPL_P_chain[row, 1], "S", IPL_S_chain[row, 1], "Z", IPL_Z_chain[row, 1])
    }
    for (col in 1:ncol(IPL_C_chain)){
      IPL_head_formulaith[, col] <- IPL_head_formulaith[, 1]
      IPL_backbone_formulaith[, col] <- IPL_backbone_formulaith[, 1]
      IPL_chain_formulaith[, col] <- IPL_chain_formulaith[, 1]
    }


    # For each sample, calculate weighted ZC from weighted elemental abundances
    weighted_IPL_formula <- IPL_C_IPL_ave # initialize vector for chain formulae
    weighted_IPL_ZC <- IPL_C_IPL_ave # initialize vector for weighted ZC
    weighted_head_formula <- IPL_C_head_ave # initialize vector for chain formulae
    weighted_head_ZC <- IPL_C_head_ave # initialize vector for weighted ZC
    weighted_backbone_formula <- IPL_C_backbone_ave # initialize vector for chain formulae
    weighted_backbone_ZC <- IPL_C_backbone_ave # initialize vector for weighted ZC
    weighted_chain_formula <- IPL_C_chain_ave # initialize vector for chain formulae
    weighted_chain_ZC <- IPL_C_chain_ave # initialize vector for weighted ZC

    default_scipen_setting <- getOption("scipen") # get default setting for scientific notation
    options(scipen=999) # turn scientific notation off temporarily
    for(sample in names(IPL_C_chain_ave)){

      sample_IPL_ave_formula <- paste0("C", signif(IPL_C_IPL_ave[sample], 3),
                                       "H", signif(IPL_H_IPL_ave[sample], 3),
                                       "O", signif(IPL_O_IPL_ave[sample], 3),
                                       "N", signif(IPL_N_IPL_ave[sample], 3),
                                       "P", signif(IPL_P_IPL_ave[sample], 3),
                                       "S", signif(IPL_S_IPL_ave[sample], 3),
                                       "Z", signif(IPL_Z_IPL_ave[sample], 3))

      sample_head_ave_formula <- paste0("C", signif(IPL_C_head_ave[sample], 3),
                                        "H", signif(IPL_H_head_ave[sample], 3),
                                        "O", signif(IPL_O_head_ave[sample], 3),
                                        "N", signif(IPL_N_head_ave[sample], 3),
                                        "P", signif(IPL_P_head_ave[sample], 3),
                                        "S", signif(IPL_S_head_ave[sample], 3),
                                        "Z", signif(IPL_Z_head_ave[sample], 3))

      sample_backbone_ave_formula <- paste0("C", signif(IPL_C_backbone_ave[sample], 3),
                                            "H", signif(IPL_H_backbone_ave[sample], 3),
                                            "O", signif(IPL_O_backbone_ave[sample], 3),
                                            "N", signif(IPL_N_backbone_ave[sample], 3),
                                            "P", signif(IPL_P_backbone_ave[sample], 3),
                                            "S", signif(IPL_S_backbone_ave[sample], 3),
                                            "Z", signif(IPL_Z_backbone_ave[sample], 3))

      sample_chain_ave_formula <- paste0("C", signif(IPL_C_chain_ave[sample], 3),
                                         "H", signif(IPL_H_chain_ave[sample], 3),
                                         "O", signif(IPL_O_chain_ave[sample], 3),
                                         "N", signif(IPL_N_chain_ave[sample], 3),
                                         "P", signif(IPL_P_chain_ave[sample], 3),
                                         "S", signif(IPL_S_chain_ave[sample], 3),
                                         "Z", signif(IPL_Z_chain_ave[sample], 3))




      weighted_IPL_formula[sample] <- sample_IPL_ave_formula

      weighted_IPL_ZC[sample] <- lipid_ZC(
        C_cc = signif(IPL_C_IPL_ave[sample], 3),
        H_cc = signif(IPL_H_IPL_ave[sample], 3),
        N_cc = signif(IPL_N_IPL_ave[sample], 3),
        O_cc = signif(IPL_O_IPL_ave[sample], 3),
        P_cc = signif(IPL_P_IPL_ave[sample], 3),
        Zplus_cc = 0,
        Zminus_cc = 0,
        S_cc = signif(IPL_S_IPL_ave[sample], 3),
        Z = signif(IPL_Z_IPL_ave[sample], 3))

      weighted_head_formula[sample] <- sample_head_ave_formula

      weighted_head_ZC[sample] <- lipid_ZC(
        C_cc = signif(IPL_C_head_ave[sample], 3),
        H_cc = signif(IPL_H_head_ave[sample], 3),
        N_cc = signif(IPL_N_head_ave[sample], 3),
        O_cc = signif(IPL_O_head_ave[sample], 3),
        P_cc = signif(IPL_P_head_ave[sample], 3),
        Zplus_cc = 0,
        Zminus_cc = 0,
        S_cc = signif(IPL_S_head_ave[sample], 3),
        Z = signif(IPL_Z_head_ave[sample], 3))

      weighted_backbone_formula[sample] <- sample_backbone_ave_formula

      weighted_backbone_ZC[sample] <- lipid_ZC(
        C_cc = signif(IPL_C_backbone_ave[sample], 3),
        H_cc = signif(IPL_H_backbone_ave[sample], 3),
        N_cc = signif(IPL_N_backbone_ave[sample], 3),
        O_cc = signif(IPL_O_backbone_ave[sample], 3),
        P_cc = signif(IPL_P_backbone_ave[sample], 3),
        Zplus_cc = 0,
        Zminus_cc = 0,
        S_cc = signif(IPL_S_backbone_ave[sample], 3),
        Z = signif(IPL_Z_backbone_ave[sample], 3))

      weighted_chain_formula[sample] <- sample_chain_ave_formula

      weighted_chain_ZC[sample] <- lipid_ZC(
        C_cc = signif(IPL_C_chain_ave[sample], 3),
        H_cc = signif(IPL_H_chain_ave[sample], 3),
        N_cc = signif(IPL_N_chain_ave[sample], 3),
        O_cc = signif(IPL_O_chain_ave[sample], 3),
        P_cc = signif(IPL_P_chain_ave[sample], 3),
        Zplus_cc = 0,
        Zminus_cc = 0,
        S_cc = signif(IPL_S_chain_ave[sample], 3),
        Z = signif(IPL_Z_chain_ave[sample], 3))
    }
    options(scipen = default_scipen_setting) # turn scientific notation back on

    IPL_nCave <- colSums(sweep(IPLareasRFmol, MARGIN = 1, nC_raw, "*")) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_nUnsatave <- colSums(sweep(IPLareasRFmol, MARGIN = 1, nUnsat_raw, "*")) / colSums(n_chains_matrix * IPLareasRFmol) ################################################## mergedIPL abundance-weighted averaging

    IPL_nPentRingave <- colSums(IPL_nPentRingith * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_nHexRingave <- colSums(IPL_nHexRingith * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_ZCave <- colSums(IPL_ZCith * IPLareasRFmol) / colSums(IPLareasRFmol)

    IPL_nCave[is.nan(IPL_nCave)] <- 0
    IPL_nUnsatave[is.nan(IPL_nUnsatave)] <- 0
    IPL_nPentRingave[is.nan(IPL_nPentRingave)] <- 0
    IPL_nHexRingave[is.nan(IPL_nHexRingave)] <- 0
    IPL_ZCave[is.nan(IPL_ZCave)] <- 0

    weighted_head_ZC[is.nan(weighted_head_ZC)] <- 0
    weighted_backbone_ZC[is.nan(weighted_backbone_ZC)] <- 0
    weighted_chain_ZC[is.nan(weighted_chain_ZC)] <- 0

    mipl[["IPL_nCave"]] <- IPL_nCave
    mipl[["IPL_nUnsatave"]] <- IPL_nUnsatave
    mipl[["IPL_nPentRingave"]] <- IPL_nPentRingave
    mipl[["IPL_nHexRingave"]] <- IPL_nHexRingave
    mipl[["IPL_ZCave"]] <- IPL_ZCave

    mipl[["IPL_formulaith"]] <- IPL_all_IPL_formulaith
    mipl[["IPL_head_formulaith"]] <- IPL_head_formulaith
    mipl[["IPL_backbone_formulaith"]] <- IPL_backbone_formulaith
    mipl[["IPL_chain_formulaith"]] <- IPL_chain_formulaith

    mipl[["weighted_IPL_formula"]] <- weighted_IPL_formula
    mipl[["weighted_IPL_ZC"]] <- weighted_IPL_ZC
    mipl[["weighted_head_formula"]] <- weighted_head_formula
    mipl[["weighted_head_ZC"]] <- weighted_head_ZC
    mipl[["weighted_backbone_formula"]] <- weighted_backbone_formula
    mipl[["weighted_backbone_ZC"]] <- weighted_backbone_ZC
    mipl[["weighted_chain_formula"]] <- weighted_chain_formula
    mipl[["weighted_chain_ZC"]] <- weighted_chain_ZC

    mipl[["IPL_C_chain_ave"]] <- IPL_C_chain_ave
    mipl[["IPL_H_chain_ave"]] <- IPL_H_chain_ave
    mipl[["IPL_O_chain_ave"]] <- IPL_O_chain_ave
    mipl[["IPL_N_chain_ave"]] <- IPL_N_chain_ave

    mipl[["C_all_IPL_matrix"]] <- C_all_IPL_matrix
    mipl[["H_all_IPL_matrix"]] <- H_all_IPL_matrix
    mipl[["N_all_IPL_matrix"]] <- N_all_IPL_matrix
    mipl[["O_all_IPL_matrix"]] <- O_all_IPL_matrix
    mipl[["S_all_IPL_matrix"]] <- S_all_IPL_matrix
    mipl[["P_all_IPL_matrix"]] <- P_all_IPL_matrix
    mipl[["Z_all_IPL_matrix"]] <- Z_all_IPL_matrix
    mipl[["C_all_head_matrix"]] <- C_all_head_matrix
    mipl[["H_all_head_matrix"]] <- H_all_head_matrix
    mipl[["N_all_head_matrix"]] <- N_all_head_matrix
    mipl[["O_all_head_matrix"]] <- O_all_head_matrix
    mipl[["S_all_head_matrix"]] <- S_all_head_matrix
    mipl[["P_all_head_matrix"]] <- P_all_head_matrix
    mipl[["Z_all_head_matrix"]] <- Z_all_head_matrix
    mipl[["C_all_backbone_matrix"]] <- C_all_backbone_matrix
    mipl[["H_all_backbone_matrix"]] <- H_all_backbone_matrix
    mipl[["N_all_backbone_matrix"]] <- N_all_backbone_matrix
    mipl[["O_all_backbone_matrix"]] <- O_all_backbone_matrix
    mipl[["S_all_backbone_matrix"]] <- S_all_backbone_matrix
    mipl[["P_all_backbone_matrix"]] <- P_all_backbone_matrix
    mipl[["Z_all_backbone_matrix"]] <- Z_all_backbone_matrix
    mipl[["C_all_chain_matrix"]] <- C_all_chain_matrix
    mipl[["H_all_chain_matrix"]] <- H_all_chain_matrix
    mipl[["N_all_chain_matrix"]] <- N_all_chain_matrix
    mipl[["O_all_chain_matrix"]] <- O_all_chain_matrix
    mipl[["S_all_chain_matrix"]] <- S_all_chain_matrix
    mipl[["P_all_chain_matrix"]] <- P_all_chain_matrix
    mipl[["Z_all_chain_matrix"]] <- Z_all_chain_matrix


    # calculate fraction of IPLs that are esters, ethers, amides, nonlinked, and hydroxylated
    y_frac_ether <- c()
    y_frac_ester <- c()
    y_frac_amide <- c()
    y_frac_nonlinkage <- c()
    y_frac_hydroxylated <- c()
    y_frac_GDGT <- c()
    y_frac_AR <- c()

    for (j in 1:ncol(mipl[["IPLareasRFmol"]])){
      y_frac_ether <- c(y_frac_ether, sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_ether_matrix"]][, j]) / sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_matrix"]][, j]))
      y_frac_ester <- c(y_frac_ester, sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_ester_matrix"]][, j]) / sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_matrix"]][, j]))
      y_frac_amide <- c(y_frac_amide, sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_amide_matrix"]][, j]) / sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_matrix"]][, j]))
      y_frac_nonlinkage <- c(y_frac_nonlinkage, sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_nonlinkage_matrix"]][, j]) / sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_matrix"]][, j]))
      y_frac_hydroxylated <- c(y_frac_hydroxylated, sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_hydroxylated_matrix"]][, j]) / sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_matrix"]][, j]))
      y_frac_GDGT <- c(y_frac_GDGT, sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_GDGT_matrix"]][, j]) / sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_matrix"]][, j]))
      y_frac_AR <- c(y_frac_AR, sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_AR_matrix"]][, j]) / sum(mipl[["IPLareasRFmol"]][, j]*mipl[["n_chains_matrix"]][, j]))
    }

    # function to remove samples that do not have peak areas for a specified IPL.
    # Used only during creation of "mergedIPL"
    rmv_smpl <- function(x){
      if(is.matrix(x) | is.data.frame(x)){
        as.matrix(x[, colSums(mipl[["IPLareasRFmol"]]) != 0, drop = FALSE])
      } else if (is.vector(x)){
        as.matrix(x[colSums(mipl[["IPLareasRFmol"]]) != 0, drop = FALSE])
      } else {
        print(paste("Error removing samples without peak areas for IPL",
          this_IPL, subcategory, subgroup, "during creation of mergedIPL"))
          # note: variable this_IPL is defined inside of IPL_process()
      }
    }

    # June 28 re-evaluation of nC calculation
    y_mean_nCith <- rmv_smpl(IPL_nCave)
    mipl[["y_mean_nCith"]] <- y_mean_nCith


    # Remove samples with no peak areas from both x, y, and mean-y values
    # whole-matrix sample removal
    y_ZCith <- rmv_smpl(mipl[["IPL_ZCith"]])
    y_nCith <- rmv_smpl(mipl[["IPL_nCith"]])
    y_nUnsatith <- rmv_smpl(mipl[["IPL_nUnsatith"]])
    y_nPentRingith <- rmv_smpl(mipl[["IPL_nPentRingith"]])
    y_nHexRingith <- rmv_smpl(mipl[["IPL_nHexRingith"]])
    wt_matrix <- rmv_smpl(IPLareasRFmol)
    IPLareas <- rmv_smpl(mipl[["IPLareas"]])
    IPL_RFith <- rmv_smpl(mipl[["IPL_RFith"]])

    IPL_formulaith <- rmv_smpl(mipl[["IPL_formulaith"]])
    IPL_head_formulaith <- rmv_smpl(mipl[["IPL_head_formulaith"]])
    IPL_backbone_formulaith <- rmv_smpl(mipl[["IPL_backbone_formulaith"]])
    IPL_chain_formulaith <- rmv_smpl(mipl[["IPL_chain_formulaith"]])


    # column-linked value (e.g. sample means) sample removal
    y_mean_ZCith <- rmv_smpl(mipl[["IPL_ZCave"]])

    weighted_IPL_formula <- rmv_smpl(mipl[["weighted_IPL_formula"]])
    weighted_IPL_ZC <- rmv_smpl(mipl[["weighted_IPL_ZC"]])
    weighted_head_formula <- rmv_smpl(mipl[["weighted_head_formula"]])
    weighted_head_ZC <- rmv_smpl(mipl[["weighted_head_ZC"]])
    weighted_backbone_formula <- rmv_smpl(mipl[["weighted_backbone_formula"]])
    weighted_backbone_ZC <- rmv_smpl(mipl[["weighted_backbone_ZC"]])
    weighted_chain_formula <- rmv_smpl(mipl[["weighted_chain_formula"]])
    weighted_chain_ZC <- rmv_smpl(mipl[["weighted_chain_ZC"]])

    IPL_C_chain_ave <- rmv_smpl(mipl[["IPL_C_chain_ave"]])
    IPL_H_chain_ave <- rmv_smpl(mipl[["IPL_H_chain_ave"]])
    IPL_O_chain_ave <- rmv_smpl(mipl[["IPL_O_chain_ave"]])
    IPL_N_chain_ave <- rmv_smpl(mipl[["IPL_N_chain_ave"]])

    C_all_IPL_matrix <- rmv_smpl(mipl[["C_all_IPL_matrix"]])
    H_all_IPL_matrix <- rmv_smpl(mipl[["H_all_IPL_matrix"]])
    N_all_IPL_matrix <- rmv_smpl(mipl[["N_all_IPL_matrix"]])
    O_all_IPL_matrix <- rmv_smpl(mipl[["O_all_IPL_matrix"]])
    S_all_IPL_matrix <- rmv_smpl(mipl[["S_all_IPL_matrix"]])
    P_all_IPL_matrix <- rmv_smpl(mipl[["P_all_IPL_matrix"]])
    Z_all_IPL_matrix <- rmv_smpl(mipl[["Z_all_IPL_matrix"]])
    C_all_head_matrix <- rmv_smpl(mipl[["C_all_head_matrix"]])
    H_all_head_matrix <- rmv_smpl(mipl[["H_all_head_matrix"]])
    N_all_head_matrix <- rmv_smpl(mipl[["N_all_head_matrix"]])
    O_all_head_matrix <- rmv_smpl(mipl[["O_all_head_matrix"]])
    S_all_head_matrix <- rmv_smpl(mipl[["S_all_head_matrix"]])
    P_all_head_matrix <- rmv_smpl(mipl[["P_all_head_matrix"]])
    Z_all_head_matrix <- rmv_smpl(mipl[["Z_all_head_matrix"]])
    C_all_backbone_matrix <- rmv_smpl(mipl[["C_all_backbone_matrix"]])
    H_all_backbone_matrix <- rmv_smpl(mipl[["H_all_backbone_matrix"]])
    N_all_backbone_matrix <- rmv_smpl(mipl[["N_all_backbone_matrix"]])
    O_all_backbone_matrix <- rmv_smpl(mipl[["O_all_backbone_matrix"]])
    S_all_backbone_matrix <- rmv_smpl(mipl[["S_all_backbone_matrix"]])
    P_all_backbone_matrix <- rmv_smpl(mipl[["P_all_backbone_matrix"]])
    Z_all_backbone_matrix <- rmv_smpl(mipl[["Z_all_backbone_matrix"]])
    C_all_chain_matrix <- rmv_smpl(mipl[["C_all_chain_matrix"]])
    H_all_chain_matrix <- rmv_smpl(mipl[["H_all_chain_matrix"]])
    N_all_chain_matrix <- rmv_smpl(mipl[["N_all_chain_matrix"]])
    O_all_chain_matrix <- rmv_smpl(mipl[["O_all_chain_matrix"]])
    S_all_chain_matrix <- rmv_smpl(mipl[["S_all_chain_matrix"]])
    P_all_chain_matrix <- rmv_smpl(mipl[["P_all_chain_matrix"]])
    Z_all_chain_matrix <- rmv_smpl(mipl[["Z_all_chain_matrix"]])

    # y_mean_nCith <- rmv_smpl(mipl[["IPL_nCave"]])
    y_mean_nUnsatith <- rmv_smpl(mipl[["IPL_nUnsatave"]])
    y_mean_nPentRingith <- rmv_smpl(mipl[["IPL_nPentRingave"]])
    y_mean_nHexRingith <- rmv_smpl(mipl[["IPL_nHexRingave"]])
    y_frac_ether <- rmv_smpl(y_frac_ether)
    y_frac_ester <- rmv_smpl(y_frac_ester)
    y_frac_amide <- rmv_smpl(y_frac_amide)
    y_frac_nonlinkage <- rmv_smpl(y_frac_nonlinkage)
    y_frac_hydroxylated <- rmv_smpl(y_frac_hydroxylated)
    y_frac_GDGT <- rmv_smpl(y_frac_GDGT)
    y_frac_AR <- rmv_smpl(y_frac_AR)

    wt_for_mean <- as.matrix(colSums(wt_matrix))

    ##### Store information about merged IPLs
    ### row-linked values

    ### column-linked values
    # abundance-weighted mean properties
    mipl[["y_mean_ZCith"]] <- y_mean_ZCith
    mipl[["IPL_formulaith"]] <- IPL_formulaith
    mipl[["IPL_head_formulaith"]] <- IPL_head_formulaith
    mipl[["IPL_backbone_formulaith"]] <- IPL_backbone_formulaith
    mipl[["IPL_chain_formulaith"]] <- IPL_chain_formulaith
    mipl[["weighted_IPL_formula"]] <- weighted_IPL_formula
    mipl[["weighted_IPL_ZC"]] <- weighted_IPL_ZC
    mipl[["weighted_head_formula"]] <- weighted_head_formula
    mipl[["weighted_head_ZC"]] <- weighted_head_ZC
    mipl[["weighted_backbone_formula"]] <- weighted_backbone_formula
    mipl[["weighted_backbone_ZC"]] <- weighted_backbone_ZC
    mipl[["weighted_chain_formula"]] <- weighted_chain_formula
    mipl[["weighted_chain_ZC"]] <- weighted_chain_ZC
    mipl[["y_mean_nCith"]] <- y_mean_nCith
    mipl[["y_mean_nUnsatith"]] <- y_mean_nUnsatith
    mipl[["y_mean_nPentRingith"]] <- y_mean_nPentRingith
    mipl[["y_mean_nHexRingith"]] <- y_mean_nHexRingith

    mipl[["IPL_C_chain_ave"]] <- IPL_C_chain_ave
    mipl[["IPL_H_chain_ave"]] <- IPL_H_chain_ave
    mipl[["IPL_O_chain_ave"]] <- IPL_O_chain_ave
    mipl[["IPL_N_chain_ave"]] <- IPL_N_chain_ave

    mipl[["C_all_IPL_matrix"]] <- C_all_IPL_matrix
    mipl[["H_all_IPL_matrix"]] <- H_all_IPL_matrix
    mipl[["N_all_IPL_matrix"]] <- N_all_IPL_matrix
    mipl[["O_all_IPL_matrix"]] <- O_all_IPL_matrix
    mipl[["S_all_IPL_matrix"]] <- S_all_IPL_matrix
    mipl[["P_all_IPL_matrix"]] <- P_all_IPL_matrix
    mipl[["Z_all_IPL_matrix"]] <- Z_all_IPL_matrix
    mipl[["C_all_head_matrix"]] <- C_all_head_matrix
    mipl[["H_all_head_matrix"]] <- H_all_head_matrix
    mipl[["N_all_head_matrix"]] <- N_all_head_matrix
    mipl[["O_all_head_matrix"]] <- O_all_head_matrix
    mipl[["S_all_head_matrix"]] <- S_all_head_matrix
    mipl[["P_all_head_matrix"]] <- P_all_head_matrix
    mipl[["Z_all_head_matrix"]] <- Z_all_head_matrix
    mipl[["C_all_backbone_matrix"]] <- C_all_backbone_matrix
    mipl[["H_all_backbone_matrix"]] <- H_all_backbone_matrix
    mipl[["N_all_backbone_matrix"]] <- N_all_backbone_matrix
    mipl[["O_all_backbone_matrix"]] <- O_all_backbone_matrix
    mipl[["S_all_backbone_matrix"]] <- S_all_backbone_matrix
    mipl[["P_all_backbone_matrix"]] <- P_all_backbone_matrix
    mipl[["Z_all_backbone_matrix"]] <- Z_all_backbone_matrix
    mipl[["C_all_chain_matrix"]] <- C_all_chain_matrix
    mipl[["H_all_chain_matrix"]] <- H_all_chain_matrix
    mipl[["N_all_chain_matrix"]] <- N_all_chain_matrix
    mipl[["O_all_chain_matrix"]] <- O_all_chain_matrix
    mipl[["S_all_chain_matrix"]] <- S_all_chain_matrix
    mipl[["P_all_chain_matrix"]] <- P_all_chain_matrix
    mipl[["Z_all_chain_matrix"]] <- Z_all_chain_matrix

    # fraction of chain linkage type
    mipl[["y_frac_ether"]] <- y_frac_ether
    mipl[["y_frac_ester"]] <- y_frac_ester
    mipl[["y_frac_amide"]] <- y_frac_amide
    mipl[["y_frac_nonlinkage"]] <- y_frac_nonlinkage
    mipl[["y_frac_hydroxylated"]] <- y_frac_hydroxylated
    mipl[["y_frac_GDGT"]] <- y_frac_GDGT
    mipl[["y_frac_AR"]] <- y_frac_AR
    # column area sums
    mipl[["wt_for_mean"]] <- wt_for_mean
    ### full matrix linked values
    mipl[["wt_matrix"]] <- wt_matrix
    mipl[["y_ZCith"]] <- y_ZCith
    mipl[["y_nCith"]] <- y_nCith
    mipl[["y_nUnsatith"]] <- y_nUnsatith
    mipl[["y_nPentRingith"]] <- y_nPentRingith
    mipl[["y_nHexRingith"]] <- y_nHexRingith
    mipl[["IPL_RFith"]] <- IPL_RFith
    mipl[["IPLareas"]] <- IPLareas
    ### list of IPLs merged in mergedIPL
    mipl[["info"]][["display"]] <- setdiff(unlist(IPL_grouping[i]), skipped_IPL)
    ### report IPL_grouping (needed by IPL_plot_func())
    mipl[["info"]][["IPL_grouping"]] <- IPL_grouping


    IPL_master[[paste("mergedIPL ", i, sep="")]] <- mipl
  }

  print(proc.time()-time)
  return(IPL_master)

}


# function to apply random noise to peak areas and response factors in a
# monte-carlo-style fashion
IPL_monte <- function(IPL_master, IPLworkbook_directory, seed,
  monte_iter = 0, peak_vary = 0.1, RF_vary = 1){

  time <- proc.time()

  set.seed(seed)

  default_scipen_setting <- getOption("scipen") # get default setting for scientific notation
  options(scipen=999) # turn scientific notation off temporarily

  # function to remove samples that do not have peak areas for a specified IPL
  # and to turn vectors into matrices readable by IPL_thermo_func
  rmv_smpl <- function(x){
    if(is.matrix(x)){
      as.matrix(x[, colSums(IPLareasRFmol) != 0, drop = FALSE])
    } else if (is.vector(x)){
      as.matrix(x[colSums(IPLareasRFmol) != 0, drop = FALSE])
    } else {
      print(paste("Error removing samples without peak areas for IPL"))
    }
  }

  # function to apply random variation to peak areas quickly
  random_area <- function(x){ x + runif(1, x*(-peak_vary), x*(peak_vary))} # summing

  IPL_master_monte <- list()
  mipl <- list()

  # Extract information from IPL_master
  IPLareas <- IPL_master[[length(IPL_master)]][["IPLareas"]]
  masses <- IPL_master[[length(IPL_master)]][["masses"]]
  IPL_RFith <- IPL_master[[length(IPL_master)]][["IPL_RFith"]]

  # IPL_names <- IPL_master[[length(IPL_master)]][["this IPL"]]
  RF_nrows <- IPL_master[[length(IPL_master)]][["IPL_RFith nrow"]]

  n_chains_matrix <- IPL_master[[length(IPL_master)]][["n_chains_matrix"]]
  n_heads_matrix <- IPL_master[[length(IPL_master)]][["n_heads_matrix"]]
  n_backbones_matrix <- IPL_master[[length(IPL_master)]][["n_backbones_matrix"]]

  C_all_IPL_matrix <- IPL_master[[length(IPL_master)]][["C_all_IPL_matrix"]]
  H_all_IPL_matrix <- IPL_master[[length(IPL_master)]][["H_all_IPL_matrix"]]
  N_all_IPL_matrix <- IPL_master[[length(IPL_master)]][["N_all_IPL_matrix"]]
  O_all_IPL_matrix <- IPL_master[[length(IPL_master)]][["O_all_IPL_matrix"]]
  S_all_IPL_matrix <- IPL_master[[length(IPL_master)]][["S_all_IPL_matrix"]]
  P_all_IPL_matrix <- IPL_master[[length(IPL_master)]][["P_all_IPL_matrix"]]
  Z_all_IPL_matrix <- IPL_master[[length(IPL_master)]][["Z_all_IPL_matrix"]]
  C_all_head_matrix <- IPL_master[[length(IPL_master)]][["C_all_head_matrix"]]
  H_all_head_matrix <- IPL_master[[length(IPL_master)]][["H_all_head_matrix"]]
  N_all_head_matrix <- IPL_master[[length(IPL_master)]][["N_all_head_matrix"]]
  O_all_head_matrix <- IPL_master[[length(IPL_master)]][["O_all_head_matrix"]]
  S_all_head_matrix <- IPL_master[[length(IPL_master)]][["S_all_head_matrix"]]
  P_all_head_matrix <- IPL_master[[length(IPL_master)]][["P_all_head_matrix"]]
  Z_all_head_matrix <- IPL_master[[length(IPL_master)]][["Z_all_head_matrix"]]
  C_all_backbone_matrix <- IPL_master[[length(IPL_master)]][["C_all_backbone_matrix"]]
  H_all_backbone_matrix <- IPL_master[[length(IPL_master)]][["H_all_backbone_matrix"]]
  N_all_backbone_matrix <- IPL_master[[length(IPL_master)]][["N_all_backbone_matrix"]]
  O_all_backbone_matrix <- IPL_master[[length(IPL_master)]][["O_all_backbone_matrix"]]
  S_all_backbone_matrix <- IPL_master[[length(IPL_master)]][["S_all_backbone_matrix"]]
  P_all_backbone_matrix <- IPL_master[[length(IPL_master)]][["P_all_backbone_matrix"]]
  Z_all_backbone_matrix <- IPL_master[[length(IPL_master)]][["Z_all_backbone_matrix"]]
  C_all_chain_matrix <- IPL_master[[length(IPL_master)]][["C_all_chain_matrix"]]
  H_all_chain_matrix <- IPL_master[[length(IPL_master)]][["H_all_chain_matrix"]]
  N_all_chain_matrix <- IPL_master[[length(IPL_master)]][["N_all_chain_matrix"]]
  O_all_chain_matrix <- IPL_master[[length(IPL_master)]][["O_all_chain_matrix"]]
  S_all_chain_matrix <- IPL_master[[length(IPL_master)]][["S_all_chain_matrix"]]
  P_all_chain_matrix <- IPL_master[[length(IPL_master)]][["P_all_chain_matrix"]]
  Z_all_chain_matrix <- IPL_master[[length(IPL_master)]][["Z_all_chain_matrix"]]

  nC_raw <- IPL_master[[length(IPL_master)]][["nC_raw"]]
  nUnsat_raw <- IPL_master[[length(IPL_master)]][["nUnsat_raw"]]
  IPL_nUnsatith <- IPL_master[[length(IPL_master)]][["IPL_nUnsatith"]]
  IPL_nPentRingith <- IPL_master[[length(IPL_master)]][["IPL_nPentRingith"]]
  IPL_nHexRingith <- IPL_master[[length(IPL_master)]][["IPL_nHexRingith"]]
  IPL_ZCith <- IPL_master[[length(IPL_master)]][["IPL_ZCith"]]

  # generate "rownums_for_RF", which keeps track of rows in mergedIPL that group each IPL to have a different RF
  rownums_for_RF <- 1
  for (i in 1:length(RF_nrows)){
    rownums_for_RF <- c(rownums_for_RF, RF_nrows[i-1] + rownums_for_RF[i-1])
  }
  RF_batch1 <- c()
  #RF_batch2 <- c()
  for (row in rownums_for_RF){
    RF_batch1 <- c(RF_batch1, IPL_RFith[row, "WAB188"]^-1) # store all RFs from batch1 in a long vector
    #RF_batch2 <- c(RF_batch2, IPL_RFith[row, "Mound OF1"]^-1) # store all RFs from batch2 in a long vector
  }

  # Begin monte carlo style iterations
  for (it in 1:monte_iter){

    IPLareas_mc <- apply(IPLareas, MARGIN = 1:2, FUN = random_area) # add random variation to peak areas

    if(length(which(IPLareas_mc < 0)) > 0){
      stop("A peak area generated by monte carlo is negative!")
    }

    IPL_RFith_mc <- IPL_RFith/IPL_RFith # initialize IPL_RFith_mc with matrix of 1s

    rand_RF_batch1 <- c()
    #rand_RF_batch2 <- c()
    for (i in 1:length(RF_batch1)){
      rand_RF_batch1[i] <- runif(1, RF_batch1[i]/RF_vary, RF_batch1[i]*RF_vary) # generate random batch1 RF within user-specified range
      #rand_RF_batch2[i] <- runif(1, RF_batch2[i]/RF_vary, RF_batch2[i]*RF_vary) # generate random batch2 RF within user-specified range
    }

    # Store random RFs in IPL_RFith_mc array
    for (i in 1:length(RF_batch1)){
      if (i < length(RF_batch1)) {
        IPL_RFith_mc[rownums_for_RF[i]:rownums_for_RF[i+1]-1, 1:6] <- rand_RF_batch1[i]
       # IPL_RFith_mc[rownums_for_RF[i]:rownums_for_RF[i+1]-1, 15:18] <- rand_RF_batch2[i]
      } else {
        IPL_RFith_mc[rownums_for_RF[i]:nrow(IPL_RFith_mc), 1:6] <- rand_RF_batch1[i]
        #IPL_RFith_mc[rownums_for_RF[i]:nrow(IPL_RFith_mc), 15:18] <- rand_RF_batch2[i]
      }
    }

    # calculate IPLareasRF (grams) from peakarea/RF
    #IPLareasRF <- IPLareas_mc / IPL_RFith_mc
    IPLareasRF <- IPLareas / IPL_RFith_mc

    # Divide grams by molar mass to calculate moles
    IPLareasRFmol <- sweep(IPLareasRF, MARGIN = 1, masses, '/')

    # calculate average elemental composition of IPLs and their components
    IPL_C_IPL_ave <- colSums(C_all_IPL_matrix * IPLareasRFmol) / colSums(IPLareasRFmol)
    IPL_C_head_ave <- colSums(C_all_head_matrix * IPLareasRFmol) / colSums(n_heads_matrix * IPLareasRFmol)
    IPL_C_backbone_ave <- colSums(C_all_backbone_matrix * IPLareasRFmol) / colSums(n_backbones_matrix * IPLareasRFmol)
    IPL_C_chain_ave <- colSums(C_all_chain_matrix * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_H_IPL_ave <- colSums(H_all_IPL_matrix * IPLareasRFmol) / colSums(IPLareasRFmol)
    IPL_H_head_ave <- colSums(H_all_head_matrix * IPLareasRFmol) / colSums(n_heads_matrix * IPLareasRFmol)
    IPL_H_backbone_ave <- colSums(H_all_backbone_matrix * IPLareasRFmol) / colSums(n_backbones_matrix * IPLareasRFmol)
    IPL_H_chain_ave <- colSums(H_all_chain_matrix * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_N_IPL_ave <- colSums(N_all_IPL_matrix * IPLareasRFmol) / colSums(IPLareasRFmol)
    IPL_N_head_ave <- colSums(N_all_head_matrix * IPLareasRFmol) / colSums(n_heads_matrix * IPLareasRFmol)
    IPL_N_backbone_ave <- colSums(N_all_backbone_matrix * IPLareasRFmol) / colSums(n_backbones_matrix * IPLareasRFmol)
    IPL_N_chain_ave <- colSums(N_all_chain_matrix * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_O_IPL_ave <- colSums(O_all_IPL_matrix * IPLareasRFmol) / colSums(IPLareasRFmol)
    IPL_O_head_ave <- colSums(O_all_head_matrix * IPLareasRFmol) / colSums(n_heads_matrix * IPLareasRFmol)
    IPL_O_backbone_ave <- colSums(O_all_backbone_matrix * IPLareasRFmol) / colSums(n_backbones_matrix * IPLareasRFmol)
    IPL_O_chain_ave <- colSums(O_all_chain_matrix * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_S_IPL_ave <- colSums(S_all_IPL_matrix * IPLareasRFmol) / colSums(IPLareasRFmol)
    IPL_S_head_ave <- colSums(S_all_head_matrix * IPLareasRFmol) / colSums(n_heads_matrix * IPLareasRFmol)
    IPL_S_backbone_ave <- colSums(S_all_backbone_matrix * IPLareasRFmol) / colSums(n_backbones_matrix * IPLareasRFmol)
    IPL_S_chain_ave <- colSums(S_all_chain_matrix * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_P_IPL_ave <- colSums(P_all_IPL_matrix * IPLareasRFmol) / colSums(IPLareasRFmol)
    IPL_P_head_ave <- colSums(P_all_head_matrix * IPLareasRFmol) / colSums(n_heads_matrix * IPLareasRFmol)
    IPL_P_backbone_ave <- colSums(P_all_backbone_matrix * IPLareasRFmol) / colSums(n_backbones_matrix * IPLareasRFmol)
    IPL_P_chain_ave <- colSums(P_all_chain_matrix * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)
    IPL_Z_IPL_ave <- colSums(Z_all_IPL_matrix * IPLareasRFmol) / colSums(IPLareasRFmol)
    IPL_Z_head_ave <- colSums(Z_all_head_matrix * IPLareasRFmol) / colSums(n_heads_matrix * IPLareasRFmol)
    IPL_Z_backbone_ave <- colSums(Z_all_backbone_matrix * IPLareasRFmol) / colSums(n_backbones_matrix * IPLareasRFmol)
    IPL_Z_chain_ave <- colSums(Z_all_chain_matrix * IPLareasRFmol) / colSums(n_chains_matrix * IPLareasRFmol)

    # For each sample, calculate weighted ZC from weighted elemental abundances
    weighted_IPL_ZC <- lipid_ZC(
      C_cc = signif(IPL_C_IPL_ave, 3),
      H_cc = signif(IPL_H_IPL_ave, 3),
      N_cc = signif(IPL_N_IPL_ave, 3),
      O_cc = signif(IPL_O_IPL_ave, 3),
      P_cc = signif(IPL_P_IPL_ave, 3),
      Zplus_cc = 0,
      Zminus_cc = 0,
      S_cc = signif(IPL_S_IPL_ave, 3),
      Z = signif(IPL_Z_IPL_ave, 3))

    weighted_head_ZC <- lipid_ZC(
      C_cc = signif(IPL_C_head_ave, 3),
      H_cc = signif(IPL_H_head_ave, 3),
      N_cc = signif(IPL_N_head_ave, 3),
      O_cc = signif(IPL_O_head_ave, 3),
      P_cc = signif(IPL_P_head_ave, 3),
      Zplus_cc = 0,
      Zminus_cc = 0,
      S_cc = signif(IPL_S_head_ave, 3),
      Z = signif(IPL_Z_head_ave, 3))

    weighted_backbone_ZC <- lipid_ZC(
      C_cc = signif(IPL_C_backbone_ave, 3),
      H_cc = signif(IPL_H_backbone_ave, 3),
      N_cc = signif(IPL_N_backbone_ave, 3),
      O_cc = signif(IPL_O_backbone_ave, 3),
      P_cc = signif(IPL_P_backbone_ave, 3),
      Zplus_cc = 0,
      Zminus_cc = 0,
      S_cc = signif(IPL_S_backbone_ave, 3),
      Z = signif(IPL_Z_backbone_ave, 3))

    weighted_chain_ZC <- lipid_ZC(
      C_cc = signif(IPL_C_chain_ave, 3),
      H_cc = signif(IPL_H_chain_ave, 3),
      N_cc = signif(IPL_N_chain_ave, 3),
      O_cc = signif(IPL_O_chain_ave, 3),
      P_cc = signif(IPL_P_chain_ave, 3),
      Zplus_cc = 0,
      Zminus_cc = 0,
      S_cc = signif(IPL_S_chain_ave, 3),
      Z = signif(IPL_Z_chain_ave, 3))
    
    
    mipl[["weighted_IPL_ZC"]] <- weighted_IPL_ZC
    mipl[["weighted_head_ZC"]] <- weighted_head_ZC
    mipl[["weighted_backbone_ZC"]] <- weighted_backbone_ZC
    mipl[["weighted_chain_ZC"]] <- weighted_chain_ZC
    
    weighted_IPL_ZC <- rmv_smpl(mipl[["weighted_IPL_ZC"]])
    weighted_head_ZC <- rmv_smpl(mipl[["weighted_head_ZC"]])
    weighted_backbone_ZC <- rmv_smpl(mipl[["weighted_backbone_ZC"]])
    weighted_chain_ZC <- rmv_smpl(mipl[["weighted_chain_ZC"]])
    
    
    rownames(weighted_IPL_ZC) <- colnames(IPLareas)
    rownames(weighted_head_ZC) <- colnames(IPLareas)
    rownames(weighted_backbone_ZC) <- colnames(IPLareas)
    rownames(weighted_chain_ZC) <- colnames(IPLareas)

    mipl[["weighted_IPL_ZC"]] <- weighted_IPL_ZC
    mipl[["weighted_head_ZC"]] <- weighted_head_ZC
    mipl[["weighted_backbone_ZC"]] <- weighted_backbone_ZC
    mipl[["weighted_chain_ZC"]] <- weighted_chain_ZC


    IPL_master_monte[[paste(it)]] <- mipl

  } # end M.C. iterations

  options(scipen = default_scipen_setting) # turn scientific notation back on
  print(proc.time()-time)
  return(IPL_master_monte)

}
