function (mode = 1:5, datadir = c(), outputdir = c(), studyname = c(), 
  f0 = 1, f1 = 0, do.report = c(2, 4, 5), configfile = c(), 
  myfun = c(), verbose = TRUE, ...) 
{
  # 检查传递给GGIR()函数的额外参数是否存在重复的参数名
  # 如果存在重复的参数名，会根据参数值是否一致给出相应的警告信息
  input = list(...)
  if (length(input) > 0) {
    if (length(input) > 1) {
      argNames = names(input)
      dupArgNames = duplicated(argNames)
      if (any(dupArgNames)) {
        for (dupi in unique(argNames[dupArgNames])) {
          dupArgValues = input[which(argNames %in% dupi)]
          if (all(dupArgValues == dupArgValues[[1]])) {
            warning(paste0("\nArgument ", dupi, " has been provided more than once. Try to avoid this."))
          }
          else {
            warning(paste0("\nArgument ", dupi, " has been provided more than once and with inconsistent values. Please fix."))
          }
        }
      }
    }
  }
  
  # 以下代码主要完成以下几个任务：修复文件路径、检查输出目录是否存在、确定要分析的文件范围（f0和f1）
  outputdir = gsub(pattern = "\\\\", replacement = "/", x = outputdir)
  datadir = gsub(pattern = "\\\\", replacement = "/", x = datadir)
  filelist = isfilelist(datadir)
  if (dir.exists(outputdir) == FALSE) 
    stop("\nDirectory specified by argument outputdir, does not exist")
  derivef0f1 = FALSE
  if (length(f0) == 0 | length(f1) == 0) {
    derivef0f1 = TRUE
  }
  else {
    if (f0 == 0 | f1 == 0) 
      derivef0f1 = TRUE
  }
  if (derivef0f1 == TRUE) {
    f0 = 1
    if (filelist == FALSE) {
      f1 <- length(dir(datadir, recursive = TRUE, ignore.case = TRUE, 
        pattern = "[.](csv|bin|Rda|wa|cw|gt3)"))
    }
    else {
      f1 = length(datadir)
    }
  }
  
  # 以下代码的主要目的是根据mode参数设置哪些分析阶段将被执行，并处理输入文件路径和输出文件路径
  dopart1 = dopart2 = dopart3 = dopart4 = dopart5 = FALSE
  if (length(which(mode == 0)) > 0) {
    dopart1 = dopart2 = dopart3 = dopart4 = dopart5 = TRUE
  }
  else {
    if (length(which(mode == 1)) > 0) 
      dopart1 = TRUE
    if (length(which(mode == 2)) > 0) 
      dopart2 = TRUE
    if (length(which(mode == 3)) > 0) 
      dopart3 = TRUE
    if (length(which(mode == 4)) > 0) 
      dopart4 = TRUE
    if (length(which(mode == 5)) > 0) 
      dopart5 = TRUE
  }
  if (length(datadir) > 0) {
    dir2fn = datadir2fnames(datadir, filelist)
    fnames = dir2fn$fnames
    fnamesfull = dir2fn$fnamesfull
  }
  if (filelist == TRUE) {
    metadatadir = paste0(outputdir, "/output_", studyname)
  }
  else {
    outputfoldername = unlist(strsplit(datadir, "/"))[length(unlist(strsplit(datadir, 
      "/")))]
    metadatadir = paste0(outputdir, "/output_", outputfoldername)
  }
  
  # 以下代码的主要目的是处理配置文件（configfile），检查其存在和格式是否正确
  configfile_csv = c()
  ex = "csv"
  if (length(configfile) > 0) {
    ex = unlist(strsplit(basename(configfile), split = "\\."))
    ex = ex[length(ex)]
  }
  if (ex == "csv") {
    if (dir.exists(metadatadir) | length(configfile) > 0) {
      if (length(configfile) > 0) {
        if (!file.exists(configfile)) {
          stop("\nDo not supply argument configfile if the configfile does not exist yet")
        }
        else {
          configfile_csv = configfile
        }
      }
      if (dir.exists(metadatadir) & length(configfile) == 
        0) {
        config_file_in_outputdir = paste0(metadatadir, 
          "/config.csv")
        if (file.exists(config_file_in_outputdir)) 
          configfile_csv = config_file_in_outputdir
      }
    }
  }
  
  # 以下代码主要执行三个任务
  # 1.从configfile_csv中提取参数并将它们分配给不同的参数列表。
  # 2.检查参数设置以确保分析可以正常进行。
  # 3.检查所需的R软件包是否已安装。
  params = extract_params(input = input, configfile_csv = configfile_csv)
  params_sleep = params$params_sleep
  params_metrics = params$params_metrics
  params_rawdata = params$params_rawdata
  params_247 = params$params_247
  params_phyact = params$params_phyact
  params_cleaning = params$params_cleaning
  params_output = params$params_output
  params_general = params$params_general
  if (params_general[["dataFormat"]] == "ukbiobank") {
    warning("\nRunnning part 3, 4, and 5 are disabled when dataFormat is ukbiobank epoch", 
      call. = FALSE)
    dopart3 = dopart4 = dopart5 = FALSE
    mode = mode[which(mode <= 2)]
  }
  if (dopart3 == TRUE & params_metrics[["do.anglez"]] == FALSE & 
    params_general[["dataFormat"]] == "raw") {
    params_metrics[["do.anglez"]] = TRUE
  }
  if (length(myfun) != 0) {
    warning("\nAre you using GGIR as online service to others? If yes, then make sure you prohibit the", 
      " user from specifying argument myfun as this poses a security risk.", 
      call. = FALSE)
    check_myfun(myfun, params_general[["windowsizes"]])
  }
  if (params_output[["visualreport"]] == TRUE & params_general[["dataFormat"]] != 
    "raw") {
    params_output[["visualreport"]] == FALSE
    warning(paste0("Turning off visualreport generation because", 
      " dataFormat is not raw."), call. = FALSE)
  }
  if (params_metrics$do.neishabouricounts == TRUE) {
    is_actilifecounts_installed = is.element("actilifecounts", 
      installed.packages()[, 1])
    if (is_actilifecounts_installed == FALSE) {
      stop("If you want to derive Neishabouricounts, please install package: actilifecounts.", 
        call. = FALSE)
    }
    else {
      if (utils::packageVersion("actilifecounts") < "1.1.0") {
        stop("Please update R package actilifecounts to version 1.1.0 or higher", 
          call. = FALSE)
      }
    }
  }
  if (params_247$cosinor == TRUE) {
    is_ActCR_installed = is.element("ActCR", installed.packages()[, 
      1])
    if (is_ActCR_installed == FALSE) {
      stop("If you want to derive circadian rhythm indicators, please install package: ActCR.", 
        call. = FALSE)
    }
  }
  checkFormat = TRUE
  if (all(dir.exists(datadir)) == TRUE) {
    rawaccfiles = dir(datadir, full.names = TRUE)[f0:f1]
  }
  else if (all(file.exists(datadir))) {
    rawaccfiles = datadir[f0:f1]
  }
  else {
    checkFormat = FALSE
  }
  if (checkFormat == TRUE) {
    is_GGIRread_installed = is.element("GGIRread", installed.packages()[, 
      1])
    is_read.gt3x_installed = is.element("read.gt3x", installed.packages()[, 
      1])
    if (is_GGIRread_installed == FALSE | is_read.gt3x_installed == 
      FALSE) {
      getExt = function(x) {
        tmp = unlist(strsplit(x, "[.]"))
        return(tmp[length(tmp)])
      }
      rawaccfiles_formats = unique(unlist(lapply(rawaccfiles, 
        FUN = getExt)))
      if (any(grepl("cwa|wav|bin", rawaccfiles_formats))) {
        if (is_GGIRread_installed == FALSE) {
          stop("If you are working with axivity, geneactiv, or genea files, please install package: GGIRread.", 
            call. = FALSE)
        }
      }
      if (any(grepl("gt3x", rawaccfiles_formats))) {
        if (is_read.gt3x_installed == FALSE) {
          stop(paste0("If you are working with actigraph files, please install package: read.gt3x.", 
            call. = FALSE))
        }
      }
    }
  }
  
  # 这段代码主要执行三个任务
  # 1.获取已安装的GGIR软件包的版本信息。
  # 2.如果verbose参数设置为TRUE，则在控制台输出版本信息、引用建议和可重现性建议。
  # 3.定义一个名为print_console_header的辅助函数，用于在控制台输出带有标题的分隔线。
  GGIRversion = "could not extract version"
  if (is.element("GGIR", installed.packages()[, 1])) {
    GGIRversion = as.character(utils::packageVersion("GGIR"))
    if (length(GGIRversion) != 1) 
      GGIRversion = sessionInfo()$otherPkgs$GGIR$Version
  }
  if (verbose == TRUE) {
    cat(paste0("\n   GGIR version: ", GGIRversion, "\n"))
    cat("\n   Do not forget to cite GGIR in your publications via a version number and\n")
    cat("   Migueles et al. 2019 JMPB. doi: 10.1123/jmpb.2018-0063. \n")
    cat("   See also: https://cran.r-project.org/package=GGIR/vignettes/GGIR.html#citing-ggir")
    cat("\n")
    cat("\n   To make your research reproducible and interpretable always report:")
    cat("\n     (1) Accelerometer brand and product name")
    cat("\n     (2) How you configured the accelerometer")
    cat("\n     (3) Study protocol and wear instructions given to the participants")
    cat("\n     (4) GGIR version")
    cat("\n     (5) How GGIR was used: Share the config.csv file or your R script.")
    cat("\n     (6) How you post-processed / cleaned GGIR output")
    cat("\n     (7) How reported outcomes relate to the specific variable names in GGIR")
  }
  print_console_header = function(headerTitle) {
    cat("\n")
    cat(paste0(rep("_", options()$width), collapse = ""))
    cat("\n", headerTitle, "\n")
  }
  
  # 这段代码负责执行GGIR的第一阶段分析，
  # 包括处理原始加速度计数据、计算基本元数据、处理时期数据以及合并多个记录。
  # 这是GGIR分析流程的第一个步骤，为后续的分析阶段提供输入数据和基本信息。
  if (dopart1 == TRUE) { 
    if (verbose == TRUE) 
      print_console_header("Part 1")
    if (!is.null(params_general[["maxRecordingInterval"]]) & 
      params_general[["overwrite"]] == TRUE) {
      basic_folder = paste0(metadatadir, "/meta/basic")
      if (dir.exists(basic_folder)) {
        basic_ms_files = dir(basic_folder, full.names = TRUE)
        if (length(basic_ms_files) > 0) {
          for (fnr in basic_ms_files) unlink(fnr, recursive = TRUE)
          rm(fnr)
        }
        rm(basic_folder, basic_ms_files)
      }
    }
    if (params_general[["dataFormat"]] == "raw") {
      g.part1(datadir = datadir, outputdir = outputdir, 
        f0 = f0, f1 = f1, studyname = studyname, myfun = myfun, 
        params_rawdata = params_rawdata, params_metrics = params_metrics, 
        params_cleaning = params_cleaning, params_general = params_general, 
        verbose = verbose)
    }
    else {
      warning(paste0("\nBe aware that you are using epoch level aggregates of raw data ", 
        "computed outside GGIR by which their reproducibility and ", 
        "transparancy is also outside the scope of GGIR. GGIR", 
        " input arguments related to raw data handling are ignored."), 
        call. = FALSE)
      convertEpochData(datadir = datadir, studyname = studyname, 
        outputdir = outputdir, params_general = params_general, 
        verbose = verbose)
    }
    if (!is.null(params_general[["maxRecordingInterval"]])) {
      appendRecords(metadatadir = metadatadir, desiredtz = params_general[["desiredtz"]], 
        idloc = params_general[["idloc"]], maxRecordingInterval = params_general[["maxRecordingInterval"]])
    }
  }
  
  
  if (dopart2 == TRUE) {
    if (verbose == TRUE) 
      print_console_header("Part 2")
    if (f1 == 0) 
      f1 = length(dir(paste0(metadatadir, "/meta/basic")))
    g.part2(datadir = datadir, metadatadir = metadatadir, 
      f0 = f0, f1 = f1, myfun = myfun, params_cleaning = params_cleaning, 
      params_247 = params_247, params_phyact = params_phyact, 
      params_output = params_output, params_general = params_general, 
      verbose = verbose)
  }
  
  
  if (dopart3 == TRUE) {
    if (verbose == TRUE) 
      print_console_header("Part 3")
    if (f1 == 0) 
      f1 = length(dir(paste0(metadatadir, "/meta/basic")))
    g.part3(metadatadir = metadatadir, f0 = f0, f1 = f1, 
      myfun = myfun, params_sleep = params_sleep, params_output = params_output, 
      params_metrics = params_metrics, params_general = params_general, 
      verbose = verbose)
  }
  if (dopart4 == TRUE) {
    if (verbose == TRUE) 
      print_console_header("Part 4")
    if (f1 == 0) 
      f1 = length(dir(paste0(metadatadir, "/meta/ms3.out")))
    g.part4(datadir = datadir, metadatadir = metadatadir, 
      f0 = f0, f1 = f1, params_sleep = params_sleep, params_metrics = params_metrics, 
      params_general = params_general, params_output = params_output, 
      params_cleaning = params_cleaning, verbose = verbose)
  }
  if (dopart5 == TRUE) {
    if (verbose == TRUE) 
      print_console_header("Part 5")
    if (f1 == 0) 
      f1 = length(dir(paste0(metadatadir, "/meta/ms3.out")))
    g.part5(datadir = datadir, metadatadir = metadatadir, 
      f0 = f0, f1 = f1, params_sleep = params_sleep, params_metrics = params_metrics, 
      params_general = params_general, params_output = params_output, 
      params_cleaning = params_cleaning, params_247 = params_247, 
      params_phyact = params_phyact, verbose = verbose)
  }
  LS = ls()
  LS = LS[which(LS %in% c("input", "txt", "derivef0f1", "dopart1", 
    "dopart2", "dopart3", "LS", "dopart4", "dopart5", "fnames", 
    "metadatadir", "ci", "config", "configfile", "filelist", 
    "outputfoldername", "numi", "logi", "conv2logical", 
    "conv2num", "SI", "params", "argNames", "dupArgNames", 
    "print_console_header", "configfile_csv", "myfun", "ex", 
    "dir2fn", "fnamesfull", "GGIRversion", "dupArgValues", 
    "verbose", "is_GGIRread_installed", "is_read.gt3x_installed", 
    "is_ActCR_installed", "is_actilifecounts_installed", 
    "rawaccfiles", "checkFormat", "getExt") == FALSE)]
  config.parameters = mget(LS)
  config.matrix = as.data.frame(createConfigFile(config.parameters, 
    GGIRversion))
  config.matrix$context[which(config.matrix$context == "")] = "not applicable"
  if (dir.exists(metadatadir)) {
    data.table::fwrite(config.matrix, file = paste0(metadatadir, 
      "/config.csv"), row.names = FALSE, sep = params_output[["sep_config"]])
  }
  else {
    if (dir.exists(datadir) == FALSE) {
      warning("\nCould not write config file because studyname or datadir are not correctly specified.")
    }
  }
  if (length(which(do.report == 2)) > 0) {
    if (verbose == TRUE) 
      print_console_header("Report part 2")
    N.files.ms2.out = length(dir(paste0(metadatadir, "/meta/ms2.out")))
    if (N.files.ms2.out > 0) {
      if (length(f0) == 0) 
        f0 = 1
      if (f1 == 0) 
        f1 = N.files.ms2.out
      if (length(params_247[["qwindow"]]) > 2 | is.character(params_247[["qwindow"]])) {
        store.long = TRUE
      }
      else {
        store.long = FALSE
      }
      g.report.part2(metadatadir = metadatadir, f0 = f0, 
        f1 = f1, maxdur = params_cleaning[["maxdur"]], 
        store.long = store.long, do.part2.pdf = params_output[["do.part2.pdf"]], 
        verbose = verbose, sep_reports = params_output[["sep_reports"]])
    }
    else {
      if (verbose == TRUE) 
        cat("\nSkipped because no milestone data available")
    }
  }
  if (length(which(do.report == 4)) > 0) {
    if (verbose == TRUE) 
      print_console_header("Report part 4")
    N.files.ms4.out = length(dir(paste0(metadatadir, "/meta/ms4.out")))
    if (N.files.ms4.out > 0) {
      if (N.files.ms4.out < f0) 
        f0 = 1
      if (N.files.ms4.out < f1) 
        f1 = N.files.ms4.out
      if (f1 == 0) 
        f1 = N.files.ms4.out
      g.report.part4(datadir = datadir, metadatadir = metadatadir, 
        f0 = f0, f1 = f1, loglocation = params_sleep[["loglocation"]], 
        storefolderstructure = params_output[["storefolderstructure"]], 
        data_cleaning_file = params_cleaning[["data_cleaning_file"]], 
        sleepwindowType = params_sleep[["sleepwindowType"]], 
        verbose = verbose, sep_reports = params_output[["sep_reports"]])
    }
    else {
      if (verbose == TRUE) 
        cat("\nSkipped because no milestone data available")
    }
  }
  if (length(which(do.report == 5)) > 0) {
    N.files.ms5.out = length(dir(paste0(metadatadir, "/meta/ms5.out")))
    if (verbose == TRUE) 
      print_console_header("Report part 5")
    if (N.files.ms5.out > 0) {
      if (N.files.ms5.out < f0) 
        f0 = 1
      if (N.files.ms5.out < f1) 
        f1 = N.files.ms5.out
      if (f1 == 0) 
        f1 = N.files.ms5.out
      g.report.part5(metadatadir = metadatadir, f0 = f0, 
        f1 = f1, loglocation = params_sleep[["loglocation"]], 
        params_cleaning = params_cleaning, week_weekend_aggregate.part5 = params_output[["week_weekend_aggregate.part5"]], 
        LUX_day_segments = params_247[["LUX_day_segments"]], 
        verbose = verbose, sep_reports = params_output[["sep_reports"]])
    }
    else {
      if (verbose == TRUE) 
        cat("\nSkipped because no milestone data available")
    }
  }
  if (params_output[["visualreport"]] == TRUE) {
    files.basic = gsub("^meta_", "", dir(paste0(metadatadir, 
      "/meta/basic")))
    files.ms3.out = dir(paste0(metadatadir, "/meta/ms3.out"))
    files.ms4.out = dir(paste0(metadatadir, "/meta/ms4.out"))
    files.available = Reduce(intersect, list(files.basic, 
      files.ms3.out, files.ms4.out))
    if (verbose == TRUE) 
      print_console_header("Generate visual reports")
    if (length(files.available) > 0) {
      g.plot5(metadatadir = metadatadir, dofirstpage = params_output[["dofirstpage"]], 
        viewingwindow = params_output[["viewingwindow"]], 
        f0 = f0, f1 = f1, overwrite = params_general[["overwrite"]], 
        metric = params_general[["acc.metric"]], desiredtz = params_general[["desiredtz"]], 
        threshold.lig = params_phyact[["threshold.lig"]], 
        threshold.mod = params_phyact[["threshold.mod"]], 
        threshold.vig = params_phyact[["threshold.vig"]], 
        visualreport_without_invalid = params_output[["visualreport_without_invalid"]], 
        includedaycrit = params_cleaning[["includedaycrit"]], 
        includenightcrit = params_cleaning[["includenightcrit"]], 
        verbose = TRUE)
    }
    else {
      if (verbose == TRUE) 
        cat("\nSkipped because no milestone data available")
    }
  }
}
