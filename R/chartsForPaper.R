
if(FALSE)
{
  library(ccInterp)
  library(dygraphs)
  library(xts)
  library(dplyr)
  library(pracma)
  library(ggplot2)
  library(data.table)
  ################################
  # start/end of data
  #

  # 9912 ringing around zero
  # 3339 just a good comparable example of bwf and cci being well matched
  #randomts <- StevesCoolRandomTS(maxFlow = 1000, obs = 600, maxNoise = 200, smoothed = TRUE, randomtimes = TRUE)


  randomSeed <- sample(1:10000, 1)
  set.seed(randomSeed)

  #set.seed(959)

  # 5698 rapid rises causes ringing at the beginning of an event.
  # 5905 rapid rises causes ringing at the beginning of an event. # good one
  randomts <- StevesCoolRandomTS(maxFlow = 200, obs = 200, maxNoise = 200, smoothed = FALSE, randomtimes = TRUE)
  #seq(randomts$Time[1], randomts$Time[1])
  # 5536 large rining before event
  randomSeed <- sample(1:10000, 1)
  set.seed(randomSeed)

  randomts <- StevesCoolRandomTS(maxFlow = 400, obs = 800, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, eventmagnitude = 0.25)


  randomts <- randomts %>% mutate(TidedSignal = Signal + Noise + TideInteraction)
  plot(randomts$Time, randomts$TidedSignal)

  randomtsHrly <- randomts %>% dplyr::select(Time, TidedSignal) %>% changeInterval(Interval = "Hourly", option = "inst")
  points(randomtsHrly, col = "red")

  bwf <- butterworthFilter(randomtsHrly)
  randomtsHrly$bwf <- bwf$Filtered
  #head(randomtsHrly)
  #head(bwf)

  cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spinterp")

  randomtsHrly$cci <- approx(cci$Date, cci$avg, randomtsHrly$Date)$y
  randomtsHrly$godin <- godinFilter(randomtsHrly$Inst)

  plot(TidedSignal ~ Time, data = randomts, type = "l", ylab = "Discharge", col = "grey")
  points(cci ~ Date, data = randomtsHrly, col = "blue", pch = 4, cex = 1)
  points(bwf ~ Date, data = randomtsHrly, col = "orange", pch = 3, cex = 1)
  points(godin ~ Date, data = randomtsHrly, col = "darkgreen", cex = 1)
  lines(Signal ~ Time, data = randomts, col = "grey15", cex = 1)

  legend("topleft", legend = c("Unfiltered", "Input Signal", "cci","bwf","godin"),
         col = c("grey", "black", "blue", "orange", "darkgreen"), lty = c(1,1,NA,NA,NA), pch = c(NA, NA, 4,3,1))

  print(randomSeed)

  cbind(
  xts(randomts$Signal, randomts$Time),
  xts(randomts$TidedSignal, randomts$Time),
  xts(randomtsHrly$bwf, randomtsHrly$Date),
  xts(randomtsHrly$cci, randomtsHrly$Date)
  ) %>% dygraph

  #####################################
  #
  # Multiple random trials
  reps <- 10000
  #bwfpvalues <- rep(0, reps)
  #ccipvalues <- rep(0, reps)
  pvalues <- list()
  rawvalues <- list()
  for(i in seq(1:reps))
  {

    # 1428 is good, 2 main picks and a dip
    # 62304 becomes worse at peak
    # 20852 is a great example, even though the cci variance is changed
    #########
    # randomSeed <- 42412 # majorly missed peak, accentuated rining at beginning of event
    #########         strong trend in balnd-altman
    randomSeed <- sample(1:1000, 1)
    randomSeed <- as.integer((as.numeric(Sys.time()) * 1000) %% 100) * randomSeed

    set.seed(randomSeed)
    print(randomSeed)

    randomts <- StevesCoolRandomTS(obs = 500, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE)

    randomts <- randomts %>% mutate(TidedSignal = Signal + Noise)
    randomtsHrly <- randomts %>% dplyr::select(Time, TidedSignal) %>% changeInterval(Interval = "Hourly", option = "inst")

    bwf <- butterworthFilter(randomtsHrly)
    randomtsHrly$bwf <- bwf$Filtered
    cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spinterp")

    randomtsHrly$cci <- approx(cci$Date, cci$avg, randomtsHrly$Date)$y
    randomtsHrly$godin <- godinFilter(randomtsHrly$Inst)

    plotdata <- randomtsHrly %>% dplyr::select(-Inst) %>%
      melt(id.vars = "Date") %>%
      na.omit %>% bind_rows(
        randomts %>% dplyr::select(c(Time, Signal, TidedSignal)) %>%
        rename("Date" = "Time") %>% melt(id.vars = "Date")
        ) %>% bind_cols(RiverTide = "No Interaction")


    # small dataframe for calculating residuals (differnce of filtered data to acutal)
    noTideInteractions <- data.frame(
      actual = approx(randomts$Time, randomts$Signal, cci$Date)$y,
      bwf = approx(bwf$Date, bwf$Filtered, cci$Date)$y,
      cci = cci$avg)

    #########################
    # With Tide Interactions
    #set.seed(1103)
    #set.seed(randomSeed)
    randomts <- randomts %>% mutate(TidedSignal = Signal + Noise + TideInteraction )
    randomtsHrly <- randomts %>% dplyr::select(Time, TidedSignal) %>% changeInterval(Interval = "Hourly", option = "inst")


    bwf <- butterworthFilter(randomtsHrly)
    randomtsHrly$bwf <- bwf$Filtered
    cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spinterp")

    randomtsHrly$cci <- approx(cci$Date, cci$avg, randomtsHrly$Date)$y
    randomtsHrly$godin <- godinFilter(randomtsHrly$Inst)


    plotdata <- bind_rows(plotdata,
      randomtsHrly %>% dplyr::select(-Inst) %>%
      melt(id.vars = "Date") %>%
      na.omit %>% bind_rows(
        randomts %>% dplyr::select(c(Time, Signal, TidedSignal)) %>%
          rename("Date" = "Time") %>% melt(id.vars = "Date")
      ) %>% bind_cols(RiverTide = "River/Tide Interaction"))



    # small dataframe for calculating residuals (differnce of filtered data to acutal)
    withTideInteractions <- data.frame(
      actual = approx(randomts$Time, randomts$Signal, cci$Date)$y,
      bwf = approx(bwf$Date, bwf$Filtered, cci$Date)$y,
      cci = cci$avg)

    ###########################
    ## Butterworth
    #############################

    res1bwf <- withTideInteractions$bwf - withTideInteractions$actual
    res2bwf <- noTideInteractions$bwf - noTideInteractions$actual
    df <- data.frame( residual_diff = res1bwf - res2bwf,
                      mean_residuals = ( res1bwf + res1bwf ) / 2 ,
                      Filter = "Butterworth")
    #######################################
    # ccInterp
    #############################
    res1cci <- withTideInteractions$cci - withTideInteractions$actual
    res2cci <- noTideInteractions$cci - noTideInteractions$actual
    df <- bind_rows(df,
                    data.frame(residual_diff = res1cci - res2cci,
                               mean_residuals = ( res1cci + res2cci ) / 2,
                               Filter = "ccInterp"))

    # dataset for levene test

    ###########################################
    # Levene test to test for significance
    # <0.05 filter result is different depending on the assumptions of synthetic tides
    library(car) # for levene

    leveneData <- bind_rows(
    withTideInteractions %>% mutate(group = "withInteractions"),
    noTideInteractions %>% mutate(group = "noInteractions"))
    modelbwf <- lm(actual ~ bwf,  data = leveneData)
    modelcci <- lm(actual ~ cci,  data = leveneData)
    levenebwf <- leveneTest(abs(residuals(modelbwf)) ~ group, data =  leveneData)
    levenecci <- leveneTest(abs(residuals(modelcci)) ~ group, data =  leveneData)

    # limit the data max to the extent of cci window
    maxsignal <- approx(randomts$Time, randomts$Signal, cci$Date)$y %>% max
    maxnoise <- approx(randomts$Time, randomts$Noise, cci$Date)$y %>% max

    pvalues[[i]] <-     data.frame(trial = i,
                                   "Butterworth" = levenebwf$`Pr(>F)`[1],
                                   "ccInterp" = levenecci$`Pr(>F)`[1]
    ) %>% melt(id.var= "trial") %>%
      mutate( TideRatio = maxsignal/maxnoise) %>%
      rename( "pvalue" = value) %>%
      mutate( d1d2Ratio = randomts$d1d2ratio[1]) %>%
      mutate( magnitude = randomts$eventmagnitude[1])

    ###################################################
    #print( paste( "BWF:", bwfpvalues[i] ) )
    #print( paste( "CCI:", ccipvalues[i] ) )

    rawvalues[[i]] <- leveneData %>% mutate( trial = i ) %>%
      mutate( TideRatio = maxsignal/maxnoise) %>%
      mutate( d1d2Ratio = randomts$d1d2ratio[1]) %>%
      mutate( magnitude = randomts$eventmagnitude[1])

    print(randomSeed)
  }


  #savedpvalues
  #pvalues <- pvaluessaved
  #rawvaluessaved <- rawvalues

  pvalues <- pvalues %>% bind_rows

  pvalues <- pvalues %>% mutate(Significant = ifelse(pvalue < 0.05, TRUE, FALSE))

  # Define refined logarithmic bins
  bin_breaks <- c(0.001, 0.0032, 0.01, 0.032, 0.1, 0.32, 1, 3.2, 10, 32, 100, Inf)
  bin_labels <- c("0.001-0.0032", "0.0032-0.01", "0.01-0.032", "0.032-0.1",
                  "0.1-0.32", "0.32-1", "1-3.2", "3.2-10", "10-32", "32-100", ">100")

  d1d2bin_breaks <- c(0.0, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99, 1.01, 1.05, 1.1, 1.25, 1.5, 2, Inf)
  d1d2_labels <- c("0-0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.9",
                   "0.9-0.95", "0.95-0.99", "0.99-1.01", "1.01-1.05",
                   "1.05-1.1", "1.1-1.25", "1.25-1.5", "1.5-2", ">2")

  # Define custom breaks based on a normal distribution centered around 0.5
  magnitude_breaks <- c(0, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1)

  # Define corresponding labels for the custom breaks
  magnitude_labels <- c("0-0.05", "0.05-0.15", "0.15-0.25", "0.25-0.35", "0.35-0.45",
                     "0.45-0.55", "0.55-0.65", "0.65-0.75", "0.75-0.85", "0.85-0.95", "0.95-1")


  pvalues <- pvalues %>%
    mutate(TideRatioBin = cut(TideRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)) %>%
    mutate(d1d2Bin = cut(d1d2Ratio , breaks = d1d2bin_breaks, labels = d1d2_labels, include.lowest = TRUE, right = FALSE)) %>%
    mutate(magnitudeBin = cut(magnitude, breaks = magnitude_breaks, labels = magnitude_labels, include.lowest = TRUE, right = FALSE))

  # Compute percentage of TRUE values per bin and filter type
  result1 <- pvalues %>%
    group_by(TideRatioBin, variable) %>%
    summarise(
      Count_True = sum(Significant),
      Total = n(),
      Percentage = (Count_True / Total) * 100,
      .groups = 'drop'
    )

  result2 <- pvalues %>%
    group_by(d1d2Bin, variable) %>%
    summarise(
      Count_True = sum(Significant),
      Total = n(),
      Percentage = (Count_True / Total) * 100,
      .groups = 'drop'
    )
  result3 <- pvalues %>%
    group_by(magnitudeBin, variable) %>%
    summarise(
      Count_True = sum(Significant),
      Total = n(),
      Percentage = (Count_True / Total) * 100,
      .groups = 'drop'
    )

  result1 %>% na.omit %>%
  # Plot results with rotated count labels
  ggplot(aes(x = TideRatioBin, y = Percentage, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = paste0("n=", Total)),
              position = position_dodge(width = 0.9),
              vjust = 0.5, hjust = -0.2,  # Adjust text positioning
              angle = 90, size = 4) +  # Rotate text 90 degrees
    #scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format as percentages
    labs(title = "Percentage of Significant Results by Event/Tide Ratio", x = "Event/Tide Ratio", y = "Percent Significant (%)", fill = "Filter") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim (c(0,100))
    coord_cartesian(clip = 'off')  # Prevent clipping of labels

  result2 %>% na.omit %>%
      # Plot results with rotated count labels
      ggplot(aes(x = d1d2Bin, y = Percentage, fill = variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label = paste0("n=", Total)),
                position = position_dodge(width = 0.9),
                vjust = 0.5, hjust = -0.2,  # Adjust text positioning
                angle = 90, size = 4) +  # Rotate text 90 degrees
      #scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format as percentages
      labs(title = "Percentage of Significant Results by d1d2 Ratio", x = "d1d2 Ratio", y = "Percent Significant (%)", fill = "Filter") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylim (c(0,100))
    coord_cartesian(clip = 'off')  # Prevent clipping of labels

  result3 %>% na.omit %>%
      # Plot results with rotated count labels
      ggplot(aes(x = magnitudeBin, y = Percentage, fill = variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label = paste0("n=", Total)),
                position = position_dodge(width = 0.9),
                vjust = 0.5, hjust = -0.2,  # Adjust text positioning
                angle = 90, size = 4) +  # Rotate text 90 degrees
      #scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format as percentages
      labs(title = "Percentage of Significant Results by EventMagnitude", x = "EventMagnitude", y = "Percent Significant (%)", fill = "Filter") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylim (c(0,100))
    coord_cartesian(clip = 'off')  # Prevent clipping of labels

  #pvaluessaved <- pvalues
  bwfpvalues <- pvalues[pvalues$variable == "Butterworth",]$pvalue
  ccipvalues <- pvalues[pvalues$variable == "ccInterp",]$pvalue

  pvalues %>% bind_rows %>% ggplot(aes(x = log(TideRatio), y = log(pvalue), colour = variable)) +
    geom_point()

    # from 1000 trials
  bwfpvalues[bwfpvalues < 0.05] %>% length
  paste( bwfpvalues[bwfpvalues < 0.05] %>% length /reps*100, "% Significant", sep = "")

  ccipvalues[ccipvalues < 0.05] %>% length
  paste( ccipvalues[ccipvalues < 0.05] %>% length /reps*100, "% Significant", sep = "")

  testingstability <- data.frame(cci = rep(0, reps), bwf = rep(0, reps))
  for(i in 1:reps)
  {
    testingstability$cci[i] <- ccipvalues[1:i][ccipvalues[1:i] < 0.05] %>% length /i * 100
    testingstability$bwf[i] <- bwfpvalues[1:i][bwfpvalues[1:i] < 0.05] %>% length /i * 100
  }
  tail(testingstability)

  testingstability$Index <- as.numeric(row.names(testingstability))
  melt(testingstability, id.vars = "Index") %>% ggplot(aes(x = Index, y = value, colour = variable)) +
    geom_line()


  rawvalues <- rawvalues %>% bind_rows

  rawvalues <- melt(rawvalues, id.vars = c('actual', 'group', 'trial', 'TideRatio','d1d2Ratio', 'magnitude' ))
  rawvalues$error <- rawvalues$value - rawvalues$actual

  rawvalues$variable

  #sample_n(rawvalues, 10000)
  anova_model <- aov(error ~ group * variable * magnitude * TideRatio * d1d2Ratio, data = rawvalues)
  summary(anova_model)

  anova_model <- aov(actual ~ value * group * variable * magnitude * TideRatio * d1d2Ratio, data = rawvalues)
  summary(anova_model)

  lm_cci_interactions <- lm(abs(error) ~ value + magnitude * log(TideRatio) * log(d1d2Ratio),
                 data = rawvalues %>% dplyr::filter(group == "withInteractions" & variable == "cci"))
  lm_cci_nointeractions <- lm(abs(error) ~ value + magnitude * log(TideRatio) * log(d1d2Ratio),
                            data = rawvalues %>% dplyr::filter(group == "noInteractions" & variable == "cci"))
  lm_bwf_interactions <- lm(abs(error) ~ value + magnitude * log(TideRatio) * log(d1d2Ratio),
                            data = rawvalues %>% dplyr::filter(group == "withInteractions" & variable == "bwf"))
  lm_bwf_nointeractions <- lm(abs(error) ~ value + magnitude * log(TideRatio) * log(d1d2Ratio),
                            data = rawvalues %>% dplyr::filter(group == "noInteractions" & variable == "bwf"))

  summary(lm_cci_interactions)
  summary(lm_bwf_interactions)
  summary(lm_cci_nointeractions)
  summary(lm_bwf_nointeractions)

  lm_cci_interactions$

  # Sample and transform
  plotdata <- rawvalues %>%
    sample_n(10000) %>%
    mutate(
      log_TideRatio = log(TideRatio),
      log_d1d2Ratio = log(d1d2Ratio)
    ) %>%
    pivot_longer(
      cols = c(magnitude, log_TideRatio, log_d1d2Ratio),
      names_to = "x_var",
      values_to = "x_value"
    )

  # Plot with facet grid and free x scales
  ggplot(plotdata, aes(x = x_value, y = abs(error), color = variable)) +
    geom_smooth(method = "loess", se = FALSE) +
    facet_grid(rows = vars(x_var), cols = vars(variable), scales = "free_x") +
    labs(
      title = "Absolute Error vs Predictor Variables by Filtering Method",
      x = "Predictor Value",
      y = "Absolute Error"
    )
  # Check assumptions: Normality of residuals
  #shapiro.test(residuals(anova_model))  # Should not be significant (p > 0.05)

  # Check assumptions: Homogeneity of variance
  leveneTest(error ~ group * variable, data = rawvalues)  # p > 0.05 means equal variance

  # Perform Tukey's HSD test for pairwise comparisons
  #tukey_result <- TukeyHSD(anova_model)
  #print(tukey_result)


  ###############################################
  # Plots for paper
  #
  # 1428 is good, 2 main picks and a dip
  # 62304 becomes worse at peak
  # randomSeed <- 20852 is a great example, even though the cci variance is changed
  #########
  randomSeed <- 875 # majorly missed peak, accentuated rining at beginning of event
  #########         strong trend in balnd-altman
  randomSeed <- sample(1:1000, 1)
  #randomSeed <- as.integer((as.numeric(Sys.time()) * 1000) %% 100) * randomSeed

  set.seed(randomSeed)
  print(randomSeed)

  par(mfrow=c(1,2))
  #set.seed(1103)
  #set.seed(randomSeed)

  # this biases towards bwf being hte worst it could possibly be
  # event/tide ratio of 2 (800 flow, 400 noise)
  # d1/d2 ratio of 1 (not significant, just middle ground)
  # event magnitude of 0.9 i.e. 90% of maximum dampening of event/tide interaction
  #randomts <- StevesCoolRandomTS(maxFlow = 800, obs = 400, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE,  d1d2ratio = c(1,1), eventmagnitude = 0.9)

  # this biases towards cci being the worst it could possibly be
  # event/tide ratio of 0.2 (80 flow, 400 noise)
  # d1/d2 ratio of 1 (not significant, just middle ground)
  # event magnitude of 0.9 i.e. 90% of maximum dampening of event/tide interaction
  randomts <- StevesCoolRandomTS(maxFlow = 80, obs = 400, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE,  d1d2ratio = c(0.95,1), eventmagnitude = 0.7)
  #randomts <- StevesCoolRandomTS(maxFlow = 400, obs = 400, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, tideInteractions = FALSE)
  #randomts <- StevesCoolRandomTS(obs = 500, smoothed = TRUE, randomtimes = TRUE, tideInteractions = FALSE)

  randomts <- randomts %>% mutate(TidedSignal = Signal + Noise)
  randomtsHrly <- randomts %>% dplyr::select(Time, TidedSignal) %>% changeInterval(Interval = "Hourly", option = "inst")

  bwf <- butterworthFilter(randomtsHrly)
  randomtsHrly$bwf <- bwf$Filtered
  cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spinterp")

  randomtsHrly$cci <- approx(cci$Date, cci$avg, randomtsHrly$Date)$y
  randomtsHrly$godin <- godinFilter(randomtsHrly$Inst)

  #randomtsHrly$Signal <- approx(randomts$Time, randomts$Signal, randomtsHrly$Date)$y

  plotdata <- randomtsHrly %>% dplyr::select(-Inst) %>%
    melt(id.vars = c("Date")) %>%
    na.omit %>% bind_rows(
      randomts %>% dplyr::select(c(Time, Signal, TidedSignal)) %>%
        rename("Date" = "Time") %>% melt(id.vars = "Date")
    ) %>% bind_cols(RiverTide = "No Interaction")




  # small dataframe for calculating residuals (differnce of filtered data to acutal)
  noTideInteractions <- data.frame(
    actual = approx(randomts$Time, randomts$Signal, cci$Date)$y,
    bwf = approx(bwf$Date, bwf$Filtered, cci$Date)$y,
    cci = cci$avg)

  #########################
  # With Tide Interactions
  #set.seed(1103)
  #set.seed(randomSeed)

  #randomts <- StevesCoolRandomTS(maxFlow = 800, obs = 400, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE)
  #randomts <- StevesCoolRandomTS(maxFlow = 400, obs = 400, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE)
  #randomts <- StevesCoolRandomTS(obs = 500, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE)

  randomts <- randomts %>% mutate(TidedSignal = Signal + Noise + TideInteraction)
  randomtsHrly <- randomts %>% dplyr::select(Time, TidedSignal) %>% changeInterval(Interval = "Hourly", option = "inst")

  bwf <- butterworthFilter(randomtsHrly)
  randomtsHrly$bwf <- bwf$Filtered
  cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spinterp")

  randomtsHrly$cci <- approx(cci$Date, cci$avg, randomtsHrly$Date)$y
  randomtsHrly$godin <- godinFilter(randomtsHrly$Inst)

  #randomtsHrly$Signal <- approx(randomts$Time, randomts$Signal, randomtsHrly$Date)$y


  plotdata <- bind_rows(plotdata,
                        randomtsHrly %>% dplyr::select(-Inst) %>%
                          melt(id.vars = c("Date")) %>%
                          na.omit %>%

                          bind_rows(
                            randomts %>% dplyr::select(c(Time, Signal, TidedSignal)) %>%
                              rename("Date" = "Time") %>% melt(id.vars = "Date")
                          ) %>% bind_cols(RiverTide = "River/Tide Interaction"))


  #View(plotdata)

  # small dataframe for calculating residuals (differnce of filtered data to acutal)
  withTideInteractions <- data.frame(
    actual = approx(randomts$Time, randomts$Signal, cci$Date)$y,
    bwf = approx(bwf$Date, bwf$Filtered, cci$Date)$y,
    cci = cci$avg)

  #xts(withTideInteractions, cci$Date) %>% dygraph
  #xts(noTideInteractions, cci$Date) %>% dygraph



  #################################################
  # simple model residual plots
  #
  bwf <- butterworthFilter(randomtsHrly)
  godin <- randomtsHrly %>% mutate(godin = godinFilter(Inst))

  cci <- ccInterpFilter(randomts %>% mutate(tided = Signal+Noise ) %>% dplyr::select(c(Time, tided)))
  #points( cci %>% dplyr::select(Date, avg) %>% changeInterval(option = "sum") , col = "blue")

  df <- data.frame(Date = cci$Date,
                   signal = approx(randomts$Time, randomts$Signal, cci$Date)$y,
                   bwf = approx(bwf$Date, bwf$Filtered, cci$Date)$y,
                   godin =  approx(godin$Date, godin$godin, cci$Date)$y,
                   cci = cci$avg)

  bwfmodel <- lm(signal ~ bwf, data = df)
  ccimodel <- lm(signal ~ cci, data = df)
  godinmodel <- lm(signal ~ godin, data = df)

  par(mfrow = c(1,1))
  plot(resid(bwfmodel) ~ fitted(bwfmodel), ylim = c(-100,100), col = "orange")
  points(resid(ccimodel) ~ fitted(ccimodel), col = "blue")
  points(resid(godinmodel) ~ fitted(godinmodel), col = "darkgreen")


  summary(bwfmodel)
  summary(ccimodel)
  summary(godinmodel)

  df


  range(randomts$Signal)
  range(randomts$Time)


  ###################
  #    # Butterworth
  #############################

  res1bwf <- withTideInteractions$bwf - withTideInteractions$actual
  res2bwf <- noTideInteractions$bwf - noTideInteractions$actual

  df <- data.frame( residual_diff = res1bwf - res2bwf,
                    mean_residuals = ( res1bwf + res1bwf ) / 2 ,
                    Filter = "Butterworth",
                    Actual = withTideInteractions$actual)
  #######################################
  # ccInterp
  #############################
  res1cci <- withTideInteractions$cci - withTideInteractions$actual
  res2cci <- noTideInteractions$cci - noTideInteractions$actual
  df <- bind_rows(df,
                  data.frame(residual_diff = res1cci - res2cci,
                             mean_residuals = ( res1cci + res2cci ) / 2,
                             Filter = "ccInterp",
                             Actual = withTideInteractions$actual))

  # dataset for levene test

  ###########################################
  # Levene test to test for significance
  # <0.05 filter result is different depending on the assumptions of synthetic tides
  library(car) # for levene

  leveneData <- bind_rows(
    withTideInteractions %>% mutate(group = "withInteractions"),
    noTideInteractions %>% mutate(group = "noInteractions"))
  modelbwf <- lm(actual ~ bwf,  data = leveneData)
  modelcci <- lm(actual ~ cci,  data = leveneData)
  levenebwf <- leveneTest(abs(residuals(modelbwf)) ~ group, data =  leveneData)
  levenecci <- leveneTest(abs(residuals(modelcci)) ~ group, data =  leveneData)
  bwfpvalues[i] <- levenebwf$`Pr(>F)`[1]
  ccipvalues[i] <- levenecci$`Pr(>F)`[1]

  modelinteractions <- lm(actual ~ value, data = melt(withTideInteractions, id.var = "actual"))
  leveneBetweenFilters <- leveneTest(abs(residuals(modelinteractions)) ~ variable, data =  melt(withTideInteractions, id.var = "actual"))
  #betweenvalues[i] <- leveneBetweenFilters$`Pr(>F)`[1]

  #betweenvaluesNoInt
  modelinteractions <- lm(actual ~ value, data = melt(noTideInteractions, id.var = "actual"))
  leveneBetweenFilters <- leveneTest(abs(residuals(modelinteractions)) ~ variable, data =  melt(noTideInteractions, id.var = "actual"))
  #betweenvaluesNoInt[i] <- leveneBetweenFilters$`Pr(>F)`[1]


  ###################################################
  print( paste( "BWF:", bwfpvalues[i] ) )
  print( paste( "CCI:", ccipvalues[i] ) )


  ##############################
  # Plotting

  ###########################
  # ggplot comparing traces
  ##########################
  color_mapping <- c(
    "TidedSignal" = "grey",
    "Signal" = "black",
    "godin" = "darkgreen",
    "cci" = "blue",
    "bwf" = "orange"
  )
  # Define line thickness
  size_mapping <- c(
    "TidedSignal" = 0.5,  # Thin
    "Signal" = 3,       # Thick
    "godin" = 1.5,        # Thin
    "cci" = 1.5,          # Thin
    "bwf" = 1.5           # Thin
  )
  plotdata$variable <- factor(plotdata$variable, levels = c("TidedSignal", "Signal", "godin", "cci", "bwf"))

  comparetraces <- plotdata %>% ggplot(aes(x = Date, y = value, colour = variable, size = variable)) +
    geom_line() +
    facet_grid(. ~ RiverTide) +
    scale_color_manual(values = color_mapping)  +
    scale_size_manual(values = size_mapping) +
    theme_minimal()
  comparetraces

  library(tidyr)
  plotdata_wide <- plotdata %>%
    pivot_wider(names_from = variable, values_from = value)
  #View(plotdata_wide)

  #######################################################
  #######################################################

  ##########################################
  #
  # Bland-Altman Plots
  #
  #######################################


  data_summary <- df %>%
    group_by(Filter) %>%
    summarise(
      mean_diff = mean(residual_diff),
      sd_diff = sd(residual_diff),
      lower_limit = mean_diff - 1.96 * sd_diff,
      upper_limit = mean_diff + 1.96 * sd_diff
    ) %>%
    mutate(p_value = case_when(
      Filter == "Butterworth" ~ bwfpvalues[i],
      Filter == "ccInterp" ~ ccipvalues[i]
    ))


  # Merge computed values back into the main data frame
  df <- merge(df, data_summary, by = "Filter")

  # Create Bland-Altman plot with facets
  blandaltman <- ggplot(df, aes(x = mean_residuals, y = residual_diff, color = Actual)) +
    geom_point() +  # Scatter plot points
    geom_hline(aes(yintercept = mean_diff), color = "red", linetype = "solid", size = 1) +  # Mean difference
    geom_hline(aes(yintercept = lower_limit), color = "red", linetype = "dashed", size = 1) +  # -1.96 SD
    geom_hline(aes(yintercept = upper_limit), color = "red", linetype = "dashed", size = 1) +  # +1.96 SD
    labs(title = "Bland-Altman Plot",
         x = "Residual mean",
         y = "Residual WithInteractions - Residual NoInteractions") +
    #theme_minimal() +
    facet_grid(. ~ Filter) +  # Facet by "Filter" (2 plots side by side)
    geom_text(data = data_summary, aes(x = min(df$mean_residuals), y = max(df$residual_diff),
                                       label = paste("p =", signif(p_value, 3))),
              hjust = 0, vjust = 1.5, fontface = "bold", inherit.aes = FALSE)

  blandaltman


library(patchwork)
comparetraces /
    blandaltman + theme_minimal()

ratiohalf /
ratiodouble


  print(randomSeed)


  ##########################################
  # End of Bland-Altman Plots
  #######################################

  #####################################
  # Response functions
  #

  randomSeed <- sample(1:1000, 1)
  randomSeed <- as.integer((as.numeric(Sys.time()) * 1000) %% 100) * randomSeed
  #randomSeed <- 1301
  set.seed(randomSeed)
  print(randomSeed)


  #########
  randomSeed <- 875 # majorly missed peak, accentuated rining at beginning of event
  #########         strong trend in balnd-altman
  randomSeed <- sample(1:1000, 1)
  #randomSeed <- as.integer((as.numeric(Sys.time()) * 1000) %% 100) * randomSeed

  set.seed(randomSeed)
  print(randomSeed)

  par(mfrow=c(1,2))
  #set.seed(1103)
  #set.seed(randomSeed)

  # this biases towards bwf being hte worst it could possibly be
  # event/tide ratio of 2 (800 flow, 400 noise)
  # d1/d2 ratio of 1 (not significant, just middle ground)
  # event magnitude of 0.9 i.e. 90% of maximum dampening of event/tide interaction
  #randomts <- StevesCoolRandomTS(maxFlow = 800, obs = 10000, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE,  d1d2ratio = c(1,1), eventmagnitude = 0.9)

  # this biases towards cci being the worst it could possibly be
  # event/tide ratio of 0.2 (80 flow, 400 noise)
  # d1/d2 ratio of 1 (not significant, just middle ground)
  # event magnitude of 0.9 i.e. 90% of maximum dampening of event/tide interaction
  randomts <- StevesCoolRandomTS(maxFlow = 80, obs = 10000, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE,  d1d2ratio = c(0.95,1), eventmagnitude = 0.7)
  #randomts <- StevesCoolRandomTS(maxFlow = 400, obs = 400, maxNoise = 400, smoothed = TRUE, randomtimes = TRUE, tideInteractions = FALSE)
  #randomts <- StevesCoolRandomTS(obs = 500, smoothed = TRUE, randomtimes = TRUE, tideInteractions = FALSE)




  #randomSeed <- 959
  #randomSeed <- 1103
  #set.seed(randomSeed)


  #randomts <- StevesCoolRandomTS(maxFlow = 1000, obs = 10000, maxNoise = 200, smoothed = TRUE, randomtimes = TRUE, tideInteractions = TRUE)

  xts(randomts$Signal+randomts$Noise+randomts$TideInteraction, randomts$Time) %>% dygraph

  plot(randomts$Time, randomts$Signal)


  randomts


  par(mfrow=c(1,2))
  plot(randomts$Time, randomts$Noise)
  plot(randomts$Time, randomts$Signal)

  randomts$Time

  # Plot the input signal
  plot(randomts$Time, randomts$Signal, type = "l", col = "blue", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Amplitude")

  randomts <- randomts %>% mutate(TidedSignal = Signal + Noise)
  plot(randomts$Time, randomts$TidedSignal)

  randomtsHrly <- randomts %>%
    dplyr::select(c(Time, TidedSignal)) %>%
    changeInterval(option = "inst", Interval = "Hourly")

  plot(TidedSignal ~ Time, data = randomts[150:200,])
  lines(randomtsHrly)


  #######################
  # apply tidal filter

  # Define a simple moving average filter
  filter_length <- 25  # Length of the moving average
  filtered_signal <- stats::filter(randomtsHrly$Inst, rep(1 / filter_length, filter_length), sides = 2)
  randomtsHrly$ma25Filter <- filtered_signal

  bwf <- butterworthFilter(randomtsHrly)
  randomtsHrly$bwf <- bwf$Filtered

  cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spinterp")
  randomtsHrly$cci <- approx(cci$Date, cci$avg, randomtsHrly$Date)$y

  # apply godin filter
  godin <- godinFilter(randomtsHrly$Inst)
  randomtsHrly$Godin <- godin

  par(mfrow = c(2,1))

  # Plot the input signal
  plot(Signal ~ Time, data = randomts, type = "l", col = "grey10", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Discharge")

  # Plot the filtered signal
  lines(bwf ~ Date, data = randomtsHrly, col = "orange", lwd = 2)
  lines(cci ~ Date, data = randomtsHrly, col = "blue", lwd = 2)
  lines(Godin ~ Date, data = randomtsHrly, col = "darkgreen", lwd = 2)
  legend("topleft", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  maxloc <- which.max(randomts$TidedSignal)
  ts_subset <- randomts[(maxloc-100):(maxloc+100),]

  # Plot the input signal
  plot(Signal ~ Time, data = ts_subset, type = "l", col = "grey10", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Discharge")

  # Plot the filtered signal
  lines(bwf ~ Date, data = randomtsHrly, col = "orange", lwd = 2)
  lines(cci ~ Date, data = randomtsHrly, col = "blue", lwd = 2)
  lines(Godin ~ Date, data = randomtsHrly, col = "darkgreen", lwd = 2)
  legend("topleft", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  # remove NA's either side
  randomtsHrly <- na.omit(randomtsHrly)

  # Spectrum of input signal
  input_spectrum <- spectrum(randomtsHrly$Inst , plot = FALSE)
  filtered_spectrum <- spectrum(randomtsHrly$bwf, plot = FALSE)

  dt <- as.numeric(median(diff(randomtsHrly$Date)))

  # Extract frequencies and amplitudes
  input_freq <- input_spectrum$freq / dt       # Adjust frequencies for sampling rate
  input_amp <- sqrt(input_spectrum$spec)       # Convert power to amplitude
  input_period <- 1/input_freq                 # convert frequency to period

  spectrum(randomtsHrly$cci)
  spectrum(randomtsHrly$bwf)
  spectrum(randomtsHrly$Inst)

  filtered_spectrum_cci <- spectrum(randomtsHrly$cci, plot = FALSE) # spectral analysis into frequency and spectral densities
  filtered_freq_cci <- filtered_spectrum_cci$freq / dt          # divide frequency by time interval
  filtered_amp_cci <- sqrt(filtered_spectrum_cci$spec)          # calculate amplitude from spectral density
  filtered_period_cci <- 1/filtered_freq_cci                    # convert from frequency to period

  filtered_spectrum_bwf <- spectrum(randomtsHrly$bwf, plot = FALSE)
  filtered_freq_bwf <- filtered_spectrum_bwf$freq / dt
  filtered_amp_bwf <- sqrt(filtered_spectrum_bwf$spec)
  filtered_period_bwf <- 1/filtered_freq_bwf

  filtered_spectrum_gdn <- spectrum(randomtsHrly$Godin, plot = FALSE)
  filtered_freq_gdn <- filtered_spectrum_gdn$freq / dt
  filtered_amp_gdn <- sqrt(filtered_spectrum_gdn$spec)
  filtered_period_gdn <- 1/filtered_freq_gdn


  par(mfrow = c(2,2))
  plot(input_period, input_amp, log = "xy", main = "input amplitude")
  plot(filtered_period_bwf, filtered_amp_bwf, log = "xy", main = "filtered amplitude (butterworth)")
  plot(filtered_period_cci, filtered_amp_cci, log = "xy", main = "filtered amplitude (cci)")
  plot(filtered_period_gdn, filtered_amp_gdn, log = "xy", main = "filtered amplitude (godin)")


  Periodograms <- bind_rows(
  data.frame (  Period = input_period, Amplitude = input_amp, Dataset = "Input" ),
  data.frame ( Period = filtered_period_bwf, Amplitude = filtered_amp_bwf, Dataset = "Butterworth" ),
  data.frame ( Period = filtered_period_cci, Amplitude = filtered_amp_cci, Dataset = "CCI" ),
  data.frame ( Period = filtered_period_gdn, Amplitude = filtered_amp_gdn, Dataset = "Godin" ))

  head(Periodograms)
  library(scales)
  Periodograms %>% ggplot(aes(x = Period, y = Amplitude)) +
    geom_point()  +
    scale_x_continuous(trans = "log2", breaks = c(3, 6, 12, 24, 48, 168, 672, 8760),
                       minor_breaks = FALSE) +  # Set appropriate minor breaks
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_minimal() +
    #annotation_logticks() +
    facet_wrap(vars(Dataset), scales = "free_y")


  # Interpolate amplitudes to match frequencies and tehn
  # Divide the filtered amplitude by the input amplitude
  response_cci <- abs(approx(x = filtered_period_cci, y = filtered_amp_cci, xout = input_period)$y) / abs(input_amp)
  response_bwf <- abs(approx(x = filtered_period_bwf, y = filtered_amp_bwf, xout = input_period)$y) / abs(input_amp)
  response_gdn <- abs(approx(x = filtered_period_gdn, y = filtered_amp_gdn, xout = input_period)$y) / abs(input_amp)

  # save results for synthetic (for later comparison to real discharge)
  input_period_synthetic <- input_period
  response_bwf_synthetic <- response_bwf
  response_cci_synthetic <- response_cci

  par(mfrow = c(1,2))

  # Plot the response function
  plot(input_period, response_bwf, type = "l", col = "orange", lwd = 1,
       main = "Response Function (log y axis)", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(0,60), log = "y")
  abline(h = 1, col = "gray", lty = 2)
  # Plot the response function
  lines(input_period, response_cci, type = "l", col = "blue", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  lines(input_period, response_gdn, type = "l", col = "darkgreen", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  legend("bottomright", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)


  head(input_period)
  tail(input_period)

  # Plot the response function
  plot(input_period, response_bwf^2, type = "l", col = "orange", lwd = 1,
       main = "Response Function", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)
  # Plot the response function
  lines(input_period, response_cci^2, type = "l", col = "blue", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  lines(input_period, response_gdn^2, type = "l", col = "darkgreen", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  legend("bottomright", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  ####################################


  plot(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$bwf, randomts$Time)$y, col = "orange", log = "xy",
       xlab = "Synthetic Input",
       ylab = "Filtered",
       main = "Comparing to synthetic input (log)"
  )
  points(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$Godin, randomts$Time)$y, col = "darkgreen")
  points(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$cci, randomts$Time)$y, col = "blue")
  legend("bottomright", legend = c("1:1", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)
  abline(0, 1)

  plot(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$bwf, randomts$Time)$y, col = "orange",
       xlab = "Synthetic Input",
       ylab = "Filtered",
       main = "Comparing to synthetic input"
  )
  points(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$Godin, randomts$Time)$y, col = "darkgreen")
  points(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$cci, randomts$Time)$y, col = "blue")
  legend("bottomright", legend = c("1:1", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  abline(0, 1)

  print(randomSeed)

  # skw to here


  library(data.table)
  library(ggplot2)
  melt(cci[3000:5000,], id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()

  maxloc <- which.max(cci$avg)
  cci_subset <- cci[(maxloc-50):(maxloc+50),]

  melt(cci_subset, id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()

  library(matrixStats)
  stdev <- rowSds(as.matrix(cci[2:(length(cci) - 1)]))

  xts(randomts$Signal, randomts$Time) %>% dygraph


  cciConfi <- data.frame(Date = cci$Date, avg = cci$avg,
                         raw = approx(randomts$Time, randomts$TidedSignal, cci$Date)$y,
                         synthetic = approx(randomts$Time, randomts$Signal, cci$Date)$y
  )
  cciConfi <- cbind(cciConfi, upper_95 = cciConfi$avg + 2*stdev)
  cciConfi <- cbind(cciConfi, lower_95 = cciConfi$avg - 2*stdev)


  melt(cciConfi[3000:5000,], id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()

  maxloc <- which.max(cciConfi$avg)
  cci_subset <- cciConfi[(maxloc-50):(maxloc+50),]

  melt(cci_subset, id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()

  maxvalue <- which.max(cciConfi$raw)
  melt(cciConfi[(maxvalue-500):(maxvalue+500),], id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()










  ##############################################################

  ######################
  ################################
  # real data

  #JRIQ <- getTSDB(client, "1120053", var = "Discharge", bucket = "tsdata3",
  #                start = as.POSIXct("2021-01-01 00:00"),
  #               stop = as.POSIXct("2024-01-01 00:00"))
  #q2024 <- RREQ[[1]] %>% dplyr::select(time, value_Discharge, QC_Discharge)
  #qPre2024 <- RREQ[[2]] %>% dplyr::select(time, value_Discharge, QC_Discharge)
  #qMerged <- mergeTS(qPre2024, q2024)
  #saveRDS(qMerged, "RREQ.rds")

  #realQ <- readRDS("data/JRIQ.rds")
  realQ <- readRDS("data/MulgraveQ.rds")

  range(realQ$time)
  range(realQ$value_Discharge)


  #randomts <- StevesCoolRandomTS()
  #randomts

  par(mfrow=c(1,2))
  plot(realQ$time, realQ$value_Discharge,
       ylab = "Johnstone Discharge (cumecs)",
       xlab = "Date",
       type = "l")


  #######################
  # apply tidal filter

  randomts <- changeInterval(data.frame(realQ$time, realQ$value_Discharge), Interval = "Hourly", option = "inst")
  names(randomts) <- c("Time", "TidedSignal")

  river_data <- ts(randomts$TidedSignal, frequency = 24*365)
  river_fft <- fft(river_data)
  amplitudes <- Mod(river_fft)  # Compute the magnitude (amplitude)
  frequencies <- seq(0, length(river_data) - 1) / length(river_data)  # Compute frequencies

  ggplot(data.frame(frequencies, amplitudes), aes(x = 1/frequencies, y = amplitudes)) +
    geom_line() +
    scale_x_log10() +  # If you want to show frequencies on a log scale
    labs(title = "Frequency Spectrum of River Data", x = "Frequency", y = "Amplitude")




  # Define a simple moving average filter
  filter_length <- 25  # Length of the moving average

  filtered_signal <- stats::filter(randomts$TidedSignal, rep(1 / filter_length, filter_length), sides = 2)
  randomts$ma25Filter <- filtered_signal

  bwf <- butterworthFilter(data.frame(randomts$Time,  randomts$TidedSignal))
  randomts$bwf <- bwf$Filtered


  cci <- ccInterpFilter(data.frame(realQ$time,  realQ$value_Discharge), type = "spinterp")
  head(cci)

  #lines(cci$Date, cci$avg, col = "blue")
  #plot(cci$Date, cci$avg)
  randomts$cci <- approx(cci$Date, cci$avg, randomts$Time)$y

  #View(randomts)

  # Plot the input signal
  plot(randomts$Time, randomts$TidedSignal, type = "l", col = "blue", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Amplitude")


  # Plot the filtered signal
  lines(randomts$Time, randomts$ma25Filter, col = "red", lwd = 2)
  lines(randomts$Time, randomts$bwf, col = "orange", lwd = 2)
  lines(cci$Date, cci$avg, col = "blue")


  godin <- godinFilter(randomts$TidedSignal)
  randomts$Godin <- godin
  randomts <- na.omit(randomts)

  # flood event peak
  randomts %>% dplyr::filter(Time > "2023-12-15" & Time < "2023-12-20") %>%
    melt(id.vars = "Time") %>%
    ggplot(aes(x= Time, y = value, colour = variable)) +
    ggtitle("Russell Event Peak") +
    geom_line()

  # Dip before rise
  randomts %>% dplyr::filter(Time > "2024-01-09" & Time < "2024-01-16") %>%
    melt(id.vars = "Time") %>%
    ggplot(aes(x= Time, y = value, colour = variable)) +
    ggtitle("Dip Before Rise") +
    geom_line()

  # Excessive ringing during low flows
  randomts %>%
    dplyr::select(-c(TidedSignal, ma25Filter)) %>%
    dplyr::filter(Time > "2022-10-02" & Time < "2022-12-16") %>%
    melt(id.vars = "Time") %>%
    ggplot(aes(x= Time, y = value, colour = variable)) +
    ggtitle("Excessive ringing during low flows") +
    #lims(y = c(0,50)) +
    geom_line()


  # Spectrum of input signal
  input_spectrum <- spectrum(randomts$TidedSignal , plot = TRUE)

  par(mfrow = c(1,1))
  plot(1/input_spectrum$freq, input_spectrum$spec, log = "xy", type = "l")



  filtered_spectrum <- spectrum(randomts$bwf, plot = TRUE)
  spectrum(randomts$cci, plot = TRUE)

  dt <- 1
  # Extract frequencies and amplitudes
  input_freq <- input_spectrum$freq / dt       # Adjust frequencies for sampling rate
  input_amp <- sqrt(input_spectrum$spec)       # Convert power to amplitude
  input_period <- 1/input_freq

  filtered_freq <- filtered_spectrum$freq / dt
  filtered_amp <- sqrt(filtered_spectrum$spec)
  filtered_period <- 1/filtered_freq


  #Interpolate amplitudes to match frequencies
  response <- approx(x = filtered_period, y = filtered_amp, xout = input_period)$y / input_amp

  # Plot the response function
  plot(input_period, response, type = "l", col = "darkgreen", lwd = 2,
       main = "Response Function", xlab = "Frequency (cycles per hour)",
       ylab = "Response (Amplitude Ratio)",
       log = "x", ylim = c(0,1))
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  ###################################

  # Spectrum of input signal
  #input_spectrum <- spectrum(randomts$TidedSignal , plot = FALSE)

  # Extract frequencies and amplitudes
  input_freq <- input_spectrum$freq / dt       # Adjust frequencies for sampling rate
  input_amp <- sqrt(input_spectrum$spec)       # Convert power to amplitude
  input_period <- 1/input_freq

  filtered_spectrum_cci <- spectrum(randomts$cci, plot = FALSE)
  filtered_freq_cci <- filtered_spectrum_cci$freq / dt
  filtered_amp_cci <- sqrt(filtered_spectrum_cci$spec)
  filtered_period_cci <- 1/filtered_freq_cci

  filtered_spectrum_bwf <- spectrum(randomts$bwf, plot = FALSE)
  filtered_freq_bwf <- filtered_spectrum_bwf$freq / dt
  filtered_amp_bwf <- sqrt(filtered_spectrum_bwf$spec)
  filtered_period_bwf <- 1/filtered_freq_bwf

  filtered_spectrum_gdn <- spectrum(randomts$Godin, plot = FALSE)
  filtered_freq_gdn <- filtered_spectrum_gdn$freq / dt
  filtered_amp_gdn <- sqrt(filtered_spectrum_gdn$spec)
  filtered_period_gdn <- 1/filtered_freq_gdn

  filtered_spectrum_ma25 <- spectrum(randomts$ma25Filter, plot = FALSE)
  filtered_freq_ma25 <- filtered_spectrum_ma25$freq / dt
  filtered_amp_ma25 <- sqrt(filtered_spectrum_ma25$spec)
  filtered_period_ma25 <- 1/filtered_freq_ma25


  #Interpolate amplitudes to match frequencies
  response_cci <- approx(x = filtered_period_cci, y = filtered_amp_cci, xout = input_period)$y / input_amp
  response_bwf <- approx(x = filtered_period_bwf, y = filtered_amp_bwf, xout = input_period)$y / input_amp
  response_gdn <- approx(x = filtered_period_gdn, y = filtered_amp_gdn, xout = input_period)$y / input_amp
  response_ma25 <- approx(x = filtered_period_ma25, y = filtered_amp_ma25, xout = input_period)$y / input_amp


  input_period_real <- input_period
  response_bwf_real <- response_bwf
  response_cci_real <- response_cci


  # Plot the response function
  plot(input_period, response_bwf, type = "l", col = "orange", lwd = 1,
       main = "Response Function (Real Discharge)", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(0,48), log = "y")

  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  # Plot the response function
  lines(input_period, response_cci, type = "l", col = "blue", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  lines(input_period, response_gdn, type = "l", col = "darkgreen", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  lines(input_period, response_ma25, type = "l", col = "darkgrey", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  # Plot the response function
  plot(input_period, response_bwf, type = "l", col = "orange", lwd = 1,
       main = "Response Function (Real Discharge)", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  # Plot the response function
  lines(input_period, response_cci, type = "l", col = "blue", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  lines(input_period, response_gdn, type = "l", col = "darkgreen", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  ####################################

  library(xts)
  library(dygraphs)

  cbind(
    xts(randomts$TidedSignal, randomts$Time),
    xts(randomts$cci, randomts$Time),
    xts(randomts$bwf, randomts$Time),
    xts(randomts$Godin, randomts$Time)
    #xts(randomts$ma25Filter, randomts$Time)
  ) %>% dygraph



  # Plot the response function
  plot(input_period_synthetic, response_bwf_synthetic, type = "l", col = "orange", lwd = 2,
       main = "Butterworth", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(0,48), log = "y")



  # Plot the response function
  points(input_period_real, response_bwf_real, type = "l", col = "blue", lwd = 2,
         main = "Butterworth", xlab = "Period (hours)",
         ylab = "Response (Amplitude Ratio)")


  # Plot the response function
  plot(input_period_synthetic, response_cci_synthetic, type = "l", col = "orange", lwd = 2,
       main = "CCI", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(0,48), log = "y")

  # Plot the response function
  points(input_period_real, response_cci_real, type = "l", col = "blue", lwd = 2,
         main = "CCI", xlab = "Period (hours)",
         ylab = "Response (Amplitude Ratio)")

  realVsSynthetic <- bind_rows(
  data.frame( Input = input_period_synthetic, Response = response_bwf_synthetic, Filter = "Butterworth", Dataset = "Synthetic" ),
  data.frame( Input = input_period_real, Response = response_bwf_real, Filter = "Butterworth", Dataset = "Real" ),
  data.frame( Input = input_period_synthetic, Response = response_cci_synthetic, Filter = "CCI", Dataset = "Synthetic" ),
  data.frame( Input = input_period_real, Response = response_cci_real, Filter = "CCI", Dataset = "Real" ))

  require(MASS) # to access Animals data sets
  require(scales) # to access break formatting functions
  #data(Animals) # load data
  realVsSynthetic %>% ggplot(aes(x = Input, y = Response, colour = Filter)) +
    geom_line()+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                   labels = trans_format("log10", math_format(10^.x)), limits = c(1,1e3)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_bw() +
    #annotation_logticks() +
    facet_grid(rows = vars(Filter), cols = vars(Dataset))
    #xlim(NA, 100)




  ############################################
  #
  # Plot the response function
  plot(input_period_synthetic, response_bwf_synthetic, type = "l", col = "orange", lwd = 2,
       main = "Butterworth", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  # Plot the response function
  points(input_period_real, response_bwf_real, type = "l", col = "blue", lwd = 2,
         main = "Butterworth", xlab = "Period (hours)",
         ylab = "Response (Amplitude Ratio)")


  # Plot the response function
  plot(input_period_synthetic, response_cci_synthetic, type = "l", col = "orange", lwd = 2,
       main = "CCI", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  # Plot the response function
  points(input_period_real, response_cci_real, type = "l", col = "blue", lwd = 2,
         main = "CCI", xlab = "Period (hours)",
         ylab = "Response (Amplitude Ratio)")



  godinsum <- changeInterval(dplyr::select(randomts, c(Time, Godin)), option = "sum")
  bwfsum <-changeInterval(dplyr::select(randomts, c(Time, bwf)), option = "sum")
  ccisum <-changeInterval(dplyr::select(randomts, c(Time, cci)), option = "sum")
  rawsum <-changeInterval(dplyr::select(randomts, c(Time, TidedSignal )), option = "sum")

  plot(rawsum$Date, bwfsum$accum)
  points(rawsum$Date,   godinsum$accum, col ="blue")
  points(rawsum$Date,   ccisum$accum, col = "red")

  ( bwfsum[nrow(bwfsum),] - rawsum[nrow(rawsum),] ) /  rawsum[nrow(rawsum),] * 100
  ( ccisum[nrow(ccisum),] - rawsum[nrow(rawsum),]  ) /  rawsum[nrow(rawsum),] * 100
  ( godinsum[nrow(godinsum),] - rawsum[nrow(rawsum),] ) /  rawsum[nrow(rawsum),] * 100
  ( rawsum[nrow(rawsum),] - rawsum[nrow(rawsum),]  ) /  rawsum[nrow(rawsum),] * 100

  "cci = 0.013%"
  "bwf = 0.014%"
  "godin = 0.016%"


  # Impulse response as per Roberts & Roberts 1978
  # To examinet he transientr esponsea, n impulseo f magnitude
  # 50 at time 200 has been applied to each of the filters. The
  # results are shown in Figure 3.

  df <- data.frame(Time = seq(Sys.time() %>% round("hour"), by = 60*60, length.out = 512),
             Q = rep(0,512))
  df[200,2] <- 50
  plot(df, ylim = c(-1,4), type = "l")

  #lines(butterworthFilter(df), col = "orange")
  bwf <- butterworthFilter(df)
  df$bwf <- bwf$Filtered
  cci <- ccInterpFilter(df)
  df$cci <- approx(cci$Date, cci$avg, df$Time)$y

  df$godin <- godinFilter(df$Q)

  df <- melt(df, id.vars = "Time")

  # impulse response
  head(df)
  df %>% ggplot(aes(x = Time, y = value, color = variable)) +
    geom_line() +
    ylim (c(-1,5))


}
