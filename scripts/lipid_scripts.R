#import best peak picking parameters from IPO
centWprams.from.file = read.csv("data/IPO/IPO_centWaveparamfits_2020-09-23T6-03-00_AM-0600.csv",colClasses = "character")

centW.min_peakwidth = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="min_peakwidth",2])
centW.max_peakwidth = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="max_peakwidth",2])
centW.ppm = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="ppm",2])
centW.mzdiff = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="mzdiff",2])
centW.snthresh = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="snthresh",2])
centW.prefilter = c(as.numeric(centWprams.from.file[centWprams.from.file[,1]=="prefilter",2]),as.numeric(centWprams.from.file[centWprams.from.file[,1]=="value_of_prefilter",2]))
centW.noise = as.numeric(centWprams.from.file[centWprams.from.file[,1]=="noise",2])
# specify some additional settings we wish to keep constant, regardless of where the parameters above were obtained
centW.fitgauss = FALSE
centW.sleep = 1
centW.mzCenterFun = c("wMean")
centW.verbose.columns = TRUE
centW.integrate = 1
centW.profparam = list(step=0.001)

#import best grouping parameters from IPO
resultRetcorGroup <- readRDS("data/IPO/resultRetcorGroup_loess.RDS")
density.bw = resultRetcorGroup$best_settings$bw
density.max = resultRetcorGroup$best_settings$max
density.minfrac = resultRetcorGroup$best_settings$minfrac
density.minsamp = resultRetcorGroup$best_settings$minsamp
density.mzwid = resultRetcorGroup$best_settings$mzwid
density.sleep = 0 


loess.missing = resultRetcorGroup$best_settings$missing
loess.extra = resultRetcorGroup$best_settings$extra
loess.span = resultRetcorGroup$best_settings$span
loess.family = resultRetcorGroup$best_settings$family

#define polarity of adducts for CAMERA
subset.polarity = "positive"

#plot up raw EICs for standards
raw_XIC_plot_standards <- function(standards, raw) {
  group_colors <- paste0(RColorBrewer::brewer.pal(8, "Spectral"))
  names(group_colors)<- c("blank", "ng_0001", "ng_001", "ng_01", "ng_1", "ng_25", "ng_5", "ng1")
  
  for (row in 1:nrow(standards)){
    XIC <- xcms::chromatogram(raw, mz = c(standards$mz_db[row]-.005, standards$mz_db[row]+.005), rt = c(standards$rt_db[row]*60 - 120, standards$rt_db[row]*60 + 120 )) 
    plot(XIC, col = group_colors[raw$sample_group], main = standards$standard[row])
  }
}

#XIC plotting for standards (processed)
XIC_plot_standards <- function(standards, data) {
  group_colors <- paste0(RColorBrewer::brewer.pal(8, "Spectral"))
  names(group_colors)<- c("blank", "ng_0001", "ng_001", "ng_01", "ng_1", "ng_25", "ng_5", "ng1")
  
  for (row in 1:nrow(standards)){
    XIC <- xcms::chromatogram(data, mz = c(standards$mz_db[row]-.005, standards$mz_db[row]+.005), rt = c(standards$rt_db[row]*60 - 60, standards$rt_db[row]*60 + 90)) 
    sample_colors <- group_colors[XIC$sample_group]
    plot(XIC, col = group_colors, lwd = 2, peakBg = sample_colors[chromPeaks(XIC)[, "sample"]], main = standards$standard[row])
  }
} 

#plot up density plots for all files  
 plot_all_ChromPeaks <- function(data, mzXMLfiles) {
      for (i in 1:length(mzXMLfiles)){
        xcms::plotChromPeaks(data, file = i)
      }
    }

 #plotting grouped peak integrations for standards
 chromdens_plot_standards <- function(standards, data) {
   group_colors <- paste0(RColorBrewer::brewer.pal(8, "Spectral"))
   names(group_colors)<- c("blank", "ng_0001", "ng_001", "ng_01", "ng_1", "ng_25", "ng_5", "ng1")

   
   pdp <- PeakDensityParam(
     sampleGroups = data$sample_group,
     bw = 2,
     minFraction = .05,
     minSamples = 1,
     binSize = .024,
     maxFeatures = 100
   )
   
   for (row in 1:nrow(standards)){
     XIC <- xcms::chromatogram(data, mz = c(standards$mz_db[row]-.005, standards$mz_db[row]+.005), rt = c(standards$rt_db[row]*60 - 60, standards$rt_db[row]*60 + 60)) 
     sample_colors <- group_colors[XIC$sample_group]
     xcms::plotChromPeakDensity(XIC, col = group_colors, param= pdp, peakBg = sample_colors[chromPeaks(XIC)[, "sample"]], main = standards$standard[row])
   }
 } 
 
# readinteger: for a given prompt, allows capture of user input as an integer; rejects non-integer input

readinteger = function(prompttext) {
  
  n = readline(prompt=prompttext)
  
  if (!grepl("^[0-9]+$", n)) {
    
    return(readinteger(prompttext))
    
  }
  
  as.integer(n)
  
}

# readyesno: for a given prompt, allows capture of user input as y or n; rejects other input

readyesno = function(prompttext) {
  
  n = readline(prompt=prompttext)
  
  if (!grepl("y|n", n)) {
    
    return(readyesno(prompttext))
    
  }
  
  as.character(n)
  
}


# verifyFileIonMode: return the ion mode of data in a particular mzXML file, by examining "polarity" attribute of each scan in the file

verifyFileIonMode = function(mzXMLfile) {
  
  rawfile = xcmsRaw(mzXMLfile) # create an xcmsraw object out of the first file
  
  # determine ion mode by examining identifier attached to scan events
  
  if (table(rawfile@polarity)["negative"]==0 & (table(rawfile@polarity)["positive"]==length(rawfile@scanindex))) { # can deduce that the file contains positive mode data
    
    filepolarity = 1 # positive
    
  } else if (table(rawfile@polarity)["positive"]==0 & (table(rawfile@polarity)["negative"]==length(rawfile@scanindex))) { # probably negative mode data
    
    filepolarity = -1 # negative
    
  } else if (table(rawfile@polarity)["positive"]>=1 & table(rawfile@polarity)["negative"]>=1) { # scans of both mode present in the file; the original .raw files weren't split by mode during initial .mzXML conversion, or something else is wrong
    
    stop("At least one file in the current dataset contains scans of more than one ion mode. Please ensure data for different ion modes have been extracted into separate files. Stopping...") # stop script if this is the case
    
  } else if (table(rawfile@polarity)["positive"]==0 & table(rawfile@polarity)["negative"]==0) {
    
    #filepolarity = 1 #all our data is positive, but should fix script
    
    stop("Can't determine ion mode of data in the first file. Check manner in which files were converted. Stopping...") # stop script if this is the case
    
  }
  
  filepolarity
  
}

# getSubsetIonMode: return the ion mode of a subset of files, using sapply of verifyFileIonMode

getSubsetIonMode = function(mzXMLfilelist) {
  
  ionmodecount = sum(sapply(mzXMLfilelist, verifyFileIonMode)) # get sum of ion mode indicators for the files in the subset
  
  if (ionmodecount==length(mzXMLfilelist)) { # can conclude that all files contain positive mode data
    
    subset.polarity = "positive"
    
  } else if (ionmodecount==-length(mzXMLfilelist)) { # can conclude that all files contain negative mode data
    
    subset.polarity = "negative"
    
  }
  
  subset.polarity
  
}

# selectXMLSubDir: allows user to choose which subset of files to process

selectXMLSubDir = function(mzXMLdirList) {
  
  print(paste0("mzXML files exist in the following directories:"))
  
  for (i in 1:length(mzXMLdirList)) {
    
    # get number of mzXML files in this directory
    numGoodFiles = length(list.files(mzXMLdirList[i], recursive = TRUE, full.names = TRUE, pattern = "*(.mzXML|.mzxml)"))
    
    if (numGoodFiles>0) { # there are .mzXML data files in this directory
      
      print(paste0(i, ". ", numGoodFiles," .mzXML files in directory '",mzXMLdirList[i],"'"))
      
    }
    
  }
  
  processDecision = readinteger("Specify which subset you'd like to process, using integer input: ")
  
  mzXMLdirList[processDecision]
  
}

# getFNmatches: returns index(es) of file names in a given file list containing the ID numbers in a match list

getFNmatches = function(filelist,IDnumlist) {
  
  unique(grep(paste(IDnumlist,collapse="|"),filelist, value=FALSE))
  
}

# genTimeStamp: generates a timestamp string based on the current system time

genTimeStamp = function () {
  
  output_DTG = format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a good format
  output_DTG = gsub(" ", "_", output_DTG) # replace any spaces
  output_DTG = gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac compatibility)
  
}

peakShape <- function(object, cor.val = 0.9, useNoise = 100) {
  require(xcms)
  ## Ensure we use adjusted retention times!
  xdata <- applyAdjustedRtime(object)
  peakmat <- chromPeaks(xdata) #extract everything
  
  res <- lapply(split.data.frame(peakmat, peakmat[, "sample"]), function(z) {
    ## Subset the data to the present file keeping only spectra that
    ## are in the peak range.
    currentSample <- filterRt(filterFile(xdata, z[1, "sample"]),
                              rt = range(z[, c("rtmin", "rtmax")]))
    ## Loading the data into memory - could be faster for some settings
    ## can however also be removed
    currentSample <- as(as(currentSample, "OnDiskMSnExp"), "MSnExp")
    corr <- numeric(nrow(z))
    for (i in seq_len(nrow(z))) {
      mzRange <- z[i, c("mzmin", "mzmax")] + c(-0.001, 0.001)
      rtRange <- z[i, c("rtmin", "rtmax")]
      suppressWarnings(
        ints <- intensity(filterMz(filterRt(currentSample, rtRange),
                                   mzRange))
      )
      ints[lengths(ints) == 0] <- 0
      ints <- as.integer(unlist(ints))
      ints <- ints[!is.na(ints)]
      ints <- ints[ints > useNoise]
      ## Only proceed if we have values left
      if (length(ints)) {
        ## Normalize/baseline subtract
        ints <- ints - min(ints)
        if (max(ints) > 0)
          ints <- ints / max(ints)
        fit <- try(nls(y ~ SSgauss(x, mu, sigma, h),
                       data.frame(x = 1:length(ints), y = ints)),
                   silent = TRUE)
        if (class(fit) == "try-error") {
          corr[i] <- 1        # Shouldn't that be 0???
        } else {
          ## calculate correlation of eics against gaussian fit
          if (sum(!is.na(ints - fitted(fit))) > 4 &&
              sum(!is.na(unique(ints))) > 4 &&
              sum(!is.na(unique(fitted(fit)))) > 4) {
            cor <- NULL
            options(show.error.messages = FALSE)
            cor <- try(cor.test(ints, fitted(fit),
                                method = "pearson",
                                use = "complete"))
            options(show.error.messages = TRUE)
            if (!is.null(cor) && cor$p.value <= 0.05)
              corr[i] <- cor$estimate
          }
        }
      }
    }                               # End for loop
    message("Peakshape evaluation: sample ",
            basename(fileNames(currentSample)), ": ", sum(corr < cor.val),
            "/", nrow(z), " peaks removed")
    z[corr >= cor.val, , drop = FALSE]
  })
  ## Need to "remember" the process history.
  ph <- processHistory(object)
  chromPeaks(object) <- do.call(rbind, res)
  ## Have to add also the process history
  object@.processHistory <- ph
  object
}

#plot_spectra_annotate
plot_spectra_annotate <- function(spectra, fragment_db, nl_db, background_db, to_include) {
  alpha <- NULL # to satisfy codetools 'no visible binding...'
  xmin <- xmax <- ymin <- ymax <- fill <- NULL # to satisfy codetools
  for (j in 1:length(spectra)) {
    #pull out masses from spectrum object
    mtc <- mz(spectra[[j]])
    #pull out intensities from spectrum object
    i <- intensity(spectra[[j]])
    #pull out feature ID
    feature_id <- spectra@elementMetadata$feature_id[j]
    #calculate 1/10 of maximum intensity in spectra
    max_10_i <- max(i)*.1
    #pull out precursor from spectrum object
    precursor <- precursorMz(spectra[[j]])
    #make dataframe of fragment masses, intensity, max intensity, precursor, and feature id
    dfr <- data.frame(i = i, mtc = mtc, feature_id = feature_id, max_10_i = max_10_i, precursor = precursor)
    #caluclate neutral loss database using precursor mass
    nl_db_j <- as.data.frame(nl_db)
    nl_db_j[, "precursor"] <- precursorMz(spectra[[j]])
    nl_db_j[, "mtc"] <-nl_db_j[, "precursor"] -nl_db_j[, "mass"]
    nl_db_j <- nl_db_j %>% filter(mtc > 0)
    #join neutral loss database with background and fragment databases
    db <- suppressMessages(bind_rows(fragment_db, nl_db_j, background_db))
    db <- db %>% select(c("mtc", "id", "type1", "type2", "diagnostic"))
    #annotate spectra dataframe with db with rolling join of .1 Da
    dfr_m <- data.table(dfr, key = "mtc")
    tm <- data.table(db, key = key(dfr_m))
    dfr_af <- tm[dfr_m, roll=.1]
    dfr_ar <- tm[dfr_m, roll=-.1]
    dfr_l <- full_join(dfr_af, dfr_ar, by = c("mtc", "i", "max_10_i", "feature_id", "precursor", "id", "type1", "type2", "diagnostic"))
    # create df of just top 50 intensities for labeling most intense fragments
    topn <- top_n(dfr_l, 100 , i)
    #remove masses within 8 Da of each other, keeping just the most intense
    topn[, "mtc_id"] <- lapply(topn[, "mtc"], function(x) c(0, cumsum(abs(diff(x) > 8))))
    topn <- topn %>% group_by(mtc_id) %>%  top_n(1, i)
    #save just intensities as vector
    topn_i <- topn$i
    background_masses <- dfr_l %>% filter(diagnostic == "background")
    background_masses <- background_masses$mtc
    #flag which masses are within 1 Da of the precursor
    dfr_l[, "diff"] <- abs(dfr_l[, "precursor"] - dfr_l[, "mtc"])
    different  <- dfr_l$diff 
    dfr_l[, "diagnostic_filt"] <- ifelse(different < 1.5, NA, dfr_l$diagnostic)
    #flag annotation from fragments that are less than 1/10 of the maximum intensity 
    i_a <- dfr_l$i
    dfr_l[, "diagnostic_filt2"] <- ifelse(i_a < max_10_i, NA, dfr_l$diagnostic_filt)
    #flag annotation of fragments that are background
    mtc_a <- dfr_l$mtc
    dfr_l[, "diagnostic_filt3"] <- ifelse(mtc_a %in% background_masses, NA, dfr_l$diagnostic_filt2)
    #select necessary columns
    dfr_l <- dfr_l %>% select(c("mtc", "i", "id", "max_10_i", "feature_id", "diagnostic_filt3"))
    #create title with precursor mass, retention time, and feature #
    #if (is.null(background_i_value) | max_i_value > (2 * background_i_value)) {
    title <- ggtitle(paste("Precursor M/Z:", 
                           round(precursorMz(spectra[[j]]), 2),
                           "; Retention time:", 
                           round(rtime(spectra[[j]]), 2),
                           "; ",
                           spectra@elementMetadata$feature_id[j]))
    #plot
    if (all(is.na(dfr_l$diagnostic_filt3))) {
      p <- ggplot(dfr_l, aes(x = mtc, xend = mtc, y = 0, yend = i)) +
        geom_segment() +
        #label masses for most abundant fragments
        geom_text(data=subset(dfr_l, i %in% topn_i & i > max_10_i*.9), aes(x= mtc, y = i + max_10_i, label = round(mtc, 2)), angle = 90) +
        scale_y_continuous(expand = expansion(mult = c(0, .2))) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(x = "M/Z", y = "Intensity")
    } else {
      p <- ggplot(dfr_l, aes(x = mtc, xend = mtc, y = 0, yend = i)) +
        geom_segment() +
        geom_segment(data=subset(dfr_l, !is.na(diagnostic_filt3)), color = "red") +
        #label masses for top 10 most abundant fragments
        geom_text(data=subset(dfr_l, i %in% topn_i & i > max_10_i*.9), aes(x= mtc, y = i + max_10_i, label = round(mtc, 2)), angle = 90) +
        #label annotated fragments
        geom_text(data=subset(dfr_l, !is.na(diagnostic_filt3) & i > max_10_i), mapping=aes(x= mtc + 5, y = i - max_10_i, label = id), color =  "red",  angle = 90)  +
        scale_y_continuous(expand = expansion(mult = c(0, .2))) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(x = "M/Z", y = "Intensity")
    }
    if (is.null(to_include) | feature_id %in% to_include)
      print(p + title)
  }
}


#ms2 annotate
ms2_annotate <-function(spectra, fragment_db, nl_db, background_db, intensity_filter) {
  #output <- data.frame(feature_id=character(), MS2_id=character()) 
  output <- list()
  #create dataframe
  for (j in 1:length(spectra)) {
    #pull out masses from spectrum object
    mtc <- mz(spectra[[j]])
    #pull out intensities from spectrum object
    i <- intensity(spectra[[j]])
    #pull out feature ID
    ms2_feature_id <- spectra@elementMetadata$feature_id[j]
    #pull out feature retention time
    ms2_rt <- rtime(spectra[[j]])
    #calculate 1/20 of maximum intensity in spectra
    max_10_i <- max(i)* intensity_filter
    #pull out precursor from spectrum object
    ms2_precursor <- precursorMz(spectra[[j]])
    #make dataframe of fragment masses, intensity, max intensity, precursor, and feature id
    dfr <- data.frame(i = i, mtc = mtc, ms2_feature_id = ms2_feature_id, max_10_i = max_10_i, ms2_precursor = ms2_precursor, ms2_rt = ms2_rt)
    #caluclate neutral loss database using precursor mass
    nl_db_j <- as.data.frame(nl_db)
    nl_db_j[, "precursor"] <- precursorMz(spectra[[j]])
    nl_db_j[, "mtc"] <-nl_db_j[, "precursor"] -nl_db_j[, "mass"]
    nl_db_j <- nl_db_j %>% filter(mtc > 0)
    #join neutral loss database with background and fragment databases
    db <- suppressMessages(bind_rows(fragment_db, nl_db_j, background_db))
    db <- db %>% select(c("mtc", "id", "type1", "type2", "diagnostic"))
    #annotate spectra dataframe with db with rolling join of .1 Da
    #dfr_m <- data.table(dfr, key = "mtc")
    #tm <- data.table(db, key = key(dfr_m))
    #dfr_af <- tm[dfr_m, roll=mass_tolerance, mult = "all"]
    # dfr_ar <- tm[dfr_m, roll=-mass_tolerance, mult = "all"]
    #dfr_a <- full_join(dfr_af, dfr_ar, by = c("mtc", "i", "max_10_i", "ms2_feature_id", "ms2_precursor", "id", "type1", "type2", "diagnostic", "ms2_rt"))
    #dfr_a <- dfr_a %>% filter(!is.na(diagnostic)) 
    dfr_a <- sqldf("select df.*, db.id, db.diagnostic
from dfr df
inner join db db on abs(db.mtc - df.mtc) < 0.1")
    #create vector of background masses 
    background_masses <- dfr_a %>% filter(diagnostic == "background")
    background_masses <- background_masses$mtc
    #flag which masses are within 1 Da of the precursor
    dfr_a[, "diff"] <- abs(dfr_a[, "ms2_precursor"] - dfr_a[, "mtc"])
    different  <- dfr_a$diff 
    dfr_a[, "diff_filt"] <- ifelse(different < 1.5, NA, "not precursor")
    #flag annotation from fragments that are less than 1/10 of the maximum intensity 
    i_a <- dfr_a$i
    dfr_a[, "intensity_filt"] <- ifelse(i_a < max_10_i, NA, "above min")
    #flag annotation of fragments that are background
    mtc_a <- dfr_a$mtc
    dfr_a[, "background_filt"] <- ifelse(mtc_a %in% background_masses, NA, "not background")
    #select necessary columns
    dfr_a <- dfr_a %>% select(c("ms2_precursor", "mtc", "ms2_rt", "ms2_feature_id", "diagnostic", "background_filt", "intensity_filt", "diff_filt"))
    #save output for binding
    output[[j]] <- dfr_a
  }
  #bind output from all iterations
  all_data = do.call(rbind, output)
  #filter out flagged intensities, background fragments, and annotations for masses close to precursor and remove non-distinct annotations
  all_data_filt <- all_data %>% mutate(id = ifelse(!is.na(background_filt) & !is.na(intensity_filt) & !is.na(diff_filt), as.character(all_data$diagnostic), "MS2 not diagnostic")) %>% select(c("ms2_feature_id", "ms2_precursor", "ms2_rt", "id", "mtc")) %>% distinct() %>% rownames_to_column("index")
  no_id <- all_data_filt %>% filter(id == "MS2 not diagnostic")
  no_id <- no_id$index
  final <- all_data_filt %>% group_by(ms2_feature_id) %>% mutate(id_1 = ifelse(length(id) > 1 & index %in% no_id, NA, as.character(all_data_filt$id))) %>% filter(!is.na(id_1)) %>% ungroup() %>% select(-c(index, id_1))
  return(final)
}
  